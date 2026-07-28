#!/usr/bin/env python3
"""
=============================================================================
01_render_videos.py  -  Carbonate-thickness videos in Winkel-Tripel projection
=============================================================================

Renders any of three videos from the 1 Myr compacted-carbonate-thickness grids:

    dm2026   DM2026 carbonate sediment thickness (batlow cmap)
    bw1991   BW1991 carbonate sediment thickness (batlow cmap)
    delta    DM2026 - BW1991 difference         (polar  cmap, symmetric range)

Common style:
    - global Winkel-Tripel projection (GMT R0/18c)
    - 200 dpi PNG frames -> ffmpeg libx264 yuv420p MP4
    - data-derived CPT range, rounded to nice steps, with `+e` end-cap arrows
      on the colorbar (background=True clamps overflow to the endpoint colour)
    - age stamp in the top-left of the figure
    - title above the map
    - colorbar below the map
    - cached per-frame PNGs so re-renders are cheap (pass --force to wipe)

The thickness videos (dm2026, bw1991) deliberately share the same CPT range
so the two MP4s can be played side by side.  That range is derived from the
combined-product max thickness.

The difference video uses a symmetric data-derived range, capped at the 99th
percentile of |DM2026 - BW1991| to keep the colour mapping useful even at
times when a single deep-sea fan inflates the absolute max.

Usage:
    cd .../steps_carbon/step8_analysis
    python 01_render_videos.py                     # all three videos
    python 01_render_videos.py --modes dm2026      # just one
    python 01_render_videos.py --modes delta --force --fps 8
    python 01_render_videos.py --cadence 5         # every 5 Myr (quick smoke test)
    python 01_render_videos.py --ages 0 60 100     # explicit time list

Runtime: ~3 s per frame, so ~8 min per video at the full 1 Myr cadence.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as cfg
# Coastlines come from paleo_continents (local-files-only Alfonso 2024
# modified Clennett-Müller plate model: both rotation files merged into a
# single pygplates RotationModel, plus the local Coastlines gpml). The
# plate-model-manager-backed paleo_coastlines module is intentionally NOT
# imported here so this script runs entirely offline against the files in
# input/Alfonso_etal_2024_modClennettMuller/.
import paleo_continents
import pygmt

# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------
PROJ = "R0/18c"
# REGION = "d" (== -Rd == -180/180/-90/90) matches the *native* longitude
# convention of the carbonate grids written by carbonate_sediment_thickness.py.
# Using REGION = "g" (-Rg == 0/360/-90/90) instead forces GMT to wrap the grid
# longitudes onto a different origin, which - combined with any pixel/gridline
# registration drift - triggers the "Longitude range too small" warning from
# grdimage. See cfg.to_gridline_nc() for the matching registration normaliser.
REGION = "d"
DPI = 200
FRAMERATE_DEFAULT = 8

CMAP_THICKNESS = "dem4"         # GMT bundled DEM-style sequential palette
CMAP_DELTA = "polar"            # diverging blue-white-red (GMT bundled)

# Fixed CPT range for the thickness videos. Hard-capped at 0..800 m so the
# colour mapping is identical across DM2026 and BW1991 (side-by-side play)
# and stable across reruns regardless of which 1 Myr frame happens to host
# the dataset-wide max.
THICKNESS_RANGE_M = (0.0, 800.0, 50.0)

DELTA_CAP_PCT = 99              # percentile to clip the symmetric delta range to

VIDEO_FRAME_DIR = cfg.OUTPUT_DIR / "video_frames"
VIDEO_FRAME_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Pre-scan helpers — set the CPT range once for the whole video so every
# frame uses the same colour mapping.
# ---------------------------------------------------------------------------
def _nice_step(span: float) -> float:
    """Pick a CPT step that yields ~20..50 intervals across the span."""
    for s in (5, 10, 20, 25, 50, 100, 200, 250, 500, 1000):
        if span / s <= 50:
            return float(s)
    return 2000.0


def compute_thickness_range(times) -> tuple[float, float, float]:
    """Scan DM2026 + BW1991 grids and return (lo, hi, step) for the
    shared thickness CPT.  Lo is always 0 (thickness is non-negative)."""
    hi = 0.0
    for t in times:
        for src in ("DM2026", "BW1991"):
            arr = cfg.load_grid(src, t).values
            m = float(np.nanmax(arr)) if np.isfinite(arr).any() else 0.0
            if m > hi:
                hi = m
    step = _nice_step(hi)
    hi = float(np.ceil(hi / step) * step)
    return 0.0, hi, step


def compute_delta_range(times) -> tuple[float, float, float]:
    """Scan the per-time difference grids and return (-cap, +cap, step)
    where cap = ceil( percentile99(|delta|) ) rounded to a nice step.
    """
    abs_max = 0.0
    for t in times:
        d = cfg.diff_grid(t).values
        d = d[np.isfinite(d)]
        if d.size == 0:
            continue
        c = float(np.percentile(np.abs(d), DELTA_CAP_PCT))
        if c > abs_max:
            abs_max = c
    step = _nice_step(2 * abs_max)
    cap = float(np.ceil(abs_max / step) * step)
    return -cap, cap, step


# ---------------------------------------------------------------------------
# Renderer
# ---------------------------------------------------------------------------
def _make_data_array(values: np.ndarray, ref: xr.DataArray) -> xr.DataArray:
    return xr.DataArray(values.astype(np.float32),
                        coords={"lat": ref["lat"].values,
                                "lon": ref["lon"].values},
                        dims=("lat", "lon"),
                        name=ref.name or "z")


def render_frame(mode: str, t: int, frame_path: Path,
                 thickness_range, delta_range):
    if mode in ("dm2026", "bw1991"):
        src = "DM2026" if mode == "dm2026" else "BW1991"
        da = cfg.load_grid(src, t)
        cpt = CMAP_THICKNESS
        lo, hi, step = thickness_range
        series = (lo, hi, step)
        cb_label = "Compacted thickness (m)"
        bg_clamp = True               # OK: thickness has no NaN/zero confusion
        # Thickness is non-negative; only the high end can overflow the CPT,
        # so add a forward (right-side) end-arrow only.
        cb_end_arrows = "+ef"
    elif mode == "delta":
        da = cfg.diff_grid(t)
        cpt = CMAP_DELTA
        lo, hi, step = delta_range
        series = (lo, hi, step)
        cb_label = "Thickness anomaly (m)"
        # background=True on a polar (diverging) cpt clobbers nan_transparent
        # in GMT 6.5 -- continent NaN cells end up painted with the bottom
        # colour of the cmap.  Leave background OFF here; the `+e` on the
        # colorbar position still gives end-arrow indicators.
        bg_clamp = False
        # Delta is symmetric and can overflow both ends, so keep both arrows.
        cb_end_arrows = "+e"
    else:
        raise ValueError(f"Unknown video mode {mode!r}")

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain",
                 FONT_ANNOT_PRIMARY="9p,Helvetica,black",
                 FONT_LABEL="10p,Helvetica,black",
                 FONT_TITLE="13p,Helvetica-Bold,black",
                 COLOR_NAN="220/220/220")

    # No map title — the time stamp (top-left) is the only annotation.
    fig.basemap(region=REGION, projection=PROJ, frame="af")

    cpt_path = "_ccd_cpt.cpt"
    pygmt.makecpt(cmap=cpt, series=list(series), continuous=True,
                  background=bg_clamp, output=cpt_path)
    # Materialise the DataArray to a temp NetCDF with *guaranteed* gridline
    # registration. Passing the xarray object directly is convenient but loses
    # control over the registration metadata that GMT reads - any pixel-
    # registered input would trigger the "Longitude range too small" warning.
    gridline_nc = cfg.to_gridline_nc(da)
    fig.grdimage(grid=str(gridline_nc), projection=PROJ, region=REGION,
                 cmap=cpt_path, nan_transparent=True)
    gridline_nc.unlink(missing_ok=True)
    # Coastlines reconstructed to age t via the local Alfonso 2024 modified
    # Clennett-Müller plate model (both rotation files merged into one
    # pygplates RotationModel, coastline geometries from the local gpml
    # under input/Alfonso_etal_2024_modClennettMuller/Coastlines/).
    # No plate-model-manager fetch.
    # The carbonate grids are themselves plotted in the present-day frame
    # (no per-cell paleo-reconstruction is applied to the raster), so the
    # overlay is a paleogeographic *reference* against which the per-pixel
    # difference field can be read. Silently skipped if pygplates isn't
    # installed or the local input files have moved.
    paleo_continents.plot_coastlines_on(fig, t,
                                        projection=PROJ, region=REGION,
                                        pen="0.4p,gray30")

    fig.colorbar(projection=PROJ, region=REGION, cmap=cpt_path,
                 frame=[f"x+l{cb_label}"],
                 position=f"JBC+w12c/0.3c+o0/1c+h{cb_end_arrows}")

    # Time stamp in the top-left corner. Y-offset reduced by 1 cm relative to
    # the previous setting (was 0.9c, now -0.1c) so the stamp sits just inside
    # the top edge of the map rather than 1 cm above it.
    fig.text(text=f"{t} Ma",
             font="22p,Helvetica-Bold,black",
             position="TL", justify="TL",
             offset="-0.15c/-0.1c",
             no_clip=True,
             region=REGION, projection=PROJ)

    fig.savefig(frame_path, dpi=DPI)
    if Path(cpt_path).exists():
        Path(cpt_path).unlink()


def stitch_video(frame_pattern: Path, out_path: Path, fps: int, n_frames: int):
    """ffmpeg PNG sequence -> MP4."""
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(fps),
        "-i", str(frame_pattern),
        "-frames:v", str(n_frames),
        "-c:v", "libx264", "-pix_fmt", "yuv420p",
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-crf", "20",
        str(out_path),
    ]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--modes", nargs="+",
                   default=["dm2026", "bw1991", "delta"],
                   choices=("dm2026", "bw1991", "delta"),
                   help="which video(s) to build (default: all three)")
    p.add_argument("--cadence", type=int, default=1,
                   help="render every N Myr (default 1)")
    p.add_argument("--ages", type=int, nargs="*",
                   help="explicit list of ages (overrides --cadence)")
    p.add_argument("--fps", type=int, default=FRAMERATE_DEFAULT)
    p.add_argument("--force", action="store_true",
                   help="wipe cached frames and re-render from scratch")
    p.add_argument("--keep-frames", action="store_true",
                   help="keep per-frame PNGs after MP4 stitch")
    p.add_argument("--skip-stitch", action="store_true",
                   help="render frames only, don't call ffmpeg")
    args = p.parse_args()

    if args.ages:
        # Render in time order (youngest -> oldest, like the paleotopo reference)
        ages = sorted(args.ages, reverse=True)
    else:
        ages = list(range(cfg.MAX_TIME_MA, cfg.MIN_TIME_MA - 1, -args.cadence))

    needs_thickness = any(m in args.modes for m in ("dm2026", "bw1991"))
    needs_delta = "delta" in args.modes

    if needs_thickness:
        # Fixed range (no dataset scan) — both thickness videos share the
        # same 0..800 m CPT so they can be played side by side.
        thickness_range = THICKNESS_RANGE_M
        print(f"\nThickness CPT range (fixed): "
              f"{thickness_range[0]:.0f} .. {thickness_range[1]:.0f} m "
              f"(step {thickness_range[2]:.0f})")
    else:
        thickness_range = (0.0, 1.0, 1.0)

    if needs_delta:
        print(f"\nScanning {len(ages)} grids to derive delta CPT range "
              f"(p{DELTA_CAP_PCT} of |dm-bw|) ...")
        delta_range = compute_delta_range(ages)
        print(f"  delta range: +/- {delta_range[1]:.0f} m "
              f"(step {delta_range[2]:.0f})")
    else:
        delta_range = (-1.0, 1.0, 1.0)

    for mode in args.modes:
        sub = VIDEO_FRAME_DIR / mode
        sub.mkdir(exist_ok=True)
        existing = sorted(sub.glob("frame_*.png"))
        if args.force and existing:
            print(f"[{mode}] --force: wiping {len(existing)} cached frame(s)")
            for f in existing:
                f.unlink()
            existing = []
        elif len(existing) > len(ages):
            print(f"[{mode}] {len(existing)} stale frames found (need {len(ages)}) — wiping")
            for f in existing:
                f.unlink()
            existing = []
        if existing and len(existing) >= len(ages):
            print(f"[{mode}] reusing {len(existing)} cached frames "
                  f"(pass --force to re-render)")

        for n, t in enumerate(ages):
            frame = sub / f"frame_{n:04d}.png"
            if frame.exists():
                continue
            print(f"[{mode}] frame {n+1}/{len(ages)} at {t} Ma -> {frame.name}")
            render_frame(mode, t, frame, thickness_range, delta_range)

        if args.skip_stitch:
            print(f"[{mode}] --skip-stitch: frames in {sub}")
            continue

        out_video = cfg.OUTPUT_DIR / f"carbonate_{mode}_{ages[0]}-{ages[-1]}Ma.mp4"
        if out_video.exists():
            out_video.unlink()
        stitch_video(sub / "frame_%04d.png", out_video, args.fps, n_frames=len(ages))
        print(f"  wrote {out_video}")

        if not args.keep_frames:
            shutil.rmtree(sub)


if __name__ == "__main__":
    main()
