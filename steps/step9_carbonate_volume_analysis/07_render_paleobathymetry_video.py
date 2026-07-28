#!/usr/bin/env python3
"""
=============================================================================
07_render_paleobathymetry_video.py  -  Paleobathymetry video (Winkel-Tripel)
=============================================================================

Renders a global Winkel-Tripel video of the Alfonso 2024 / pyBacktrack merged
paleobathymetry grids (one frame per 1 Myr), with reconstructed continental
polygon outlines overlaid via the local Alfonso-modified Clennett-Müller
plate model (`paleo_continents.py`).

Style mirrors `01_render_videos.py`:
    - global Winkel-Tripel projection (GMT R0/18c)
    - 200 dpi PNG frames -> ffmpeg libx264 yuv420p MP4
    - fixed CPT range (no per-frame scan; same range for every frame)
    - time stamp in the top-left corner (no map title)
    - colorbar below the map
    - cached per-frame PNGs (pass --force to wipe)
    - gridline registration forced on every grid before grdimage to avoid
      the "Longitude range too small" warning

Default colour mapping: GMT `geo` palette, range -7000 .. 0 m, step 500 m,
with a single low-end arrow on the colorbar (deepest bathymetry overflow).
Override with --cmap / --vmin / --vmax / --step if a different mapping is
wanted.

Usage:
    cd .../steps_carbon/step8_analysis
    python 07_render_paleobathymetry_video.py
    python 07_render_paleobathymetry_video.py --cadence 5
    python 07_render_paleobathymetry_video.py --ages 0 60 100
    python 07_render_paleobathymetry_video.py --force --fps 8

Runtime: ~3 s per frame, ~8 min at the full 1 Myr cadence.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as cfg
# paleo_continents now exposes BOTH continental-polygon and coastline
# reconstruction (local-files-only, both rotation files merged into one
# RotationModel). The plate-model-manager-backed paleo_coastlines module
# is intentionally NOT imported here so this script runs entirely offline
# against the files in input/Alfonso_etal_2024_modClennettMuller/.
import paleo_continents
import pygmt

# ---------------------------------------------------------------------------
# Style constants (mirroring 01_render_videos.py)
# ---------------------------------------------------------------------------
PROJ = "R0/18c"
REGION = "d"                    # -180/180/-90/90, matches paleobath lon convention
DPI = 200
FRAMERATE_DEFAULT = 8

# Default paleobathymetry colour mapping. `geo` is GMT's standard
# topo/bathy palette: deep blues at low values, light cyan near sea level.
# Restricting the upper bound to 0 m keeps the colour ramp focused on the
# ocean part of the signal. The lower bound is fixed at -6000 m: anything
# deeper (trenches) is clamped to the -6000 m colour by `background=True`
# on makecpt so the dynamic range of the ramp isn't wasted on a handful of
# outlier cells. The `+eb` arrow on the colorbar tells the reader that
# values beyond -6000 m exist and have been clamped.
DEFAULT_CMAP = "geo"
DEFAULT_VMIN, DEFAULT_VMAX, DEFAULT_STEP = -6000.0, 0.0, 500.0

# Continental polygon outline pen. Red lines, no fill (fig.plot draws the
# (lon, lat) sequence as a polyline so polygons appear as outlines only -
# no `fill=` argument is passed downstream). Slightly thicker than the
# coastline pen so the two overlays read as a clear hierarchy: red
# continental crust extent first, green present-day coastline on top.
CONT_POLY_PEN = "0.6p,red"

# Reconstructed paleo-coastline pen. Green so it sits visibly on top of
# the red continental polygons.
COASTLINE_PEN = "0.4p,green"

VIDEO_FRAME_DIR = cfg.OUTPUT_DIR / "video_frames"
VIDEO_FRAME_DIR.mkdir(parents=True, exist_ok=True)
MODE = "paleobath"


# ---------------------------------------------------------------------------
# Renderer
# ---------------------------------------------------------------------------
def render_frame(t: int, frame_path: Path, cmap: str,
                 vmin: float, vmax: float, step: float):
    da = cfg.load_paleobath(t)

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain",
                 FONT_ANNOT_PRIMARY="9p,Helvetica,black",
                 FONT_LABEL="10p,Helvetica,black",
                 FONT_TITLE="13p,Helvetica-Bold,black",
                 COLOR_NAN="220/220/220")

    # No map title - the time stamp (top-left) is the only annotation.
    fig.basemap(region=REGION, projection=PROJ, frame="af")

    cpt_path = "_paleobath_cpt.cpt"
    # background=True clamps overflow values to the CPT endpoint colours,
    # which is what we want for paleobathymetry (deeper-than-vmin trenches
    # clamp to the deepest blue rather than rendering as NaN/grey).
    pygmt.makecpt(cmap=cmap, series=[vmin, vmax, step],
                  continuous=True, background=True, output=cpt_path)

    # Force gridline registration to silence the "Longitude range too small"
    # warning regardless of how the upstream paleobathymetry was registered.
    gridline_nc = cfg.to_gridline_nc(da)
    fig.grdimage(grid=str(gridline_nc),
                 projection=PROJ, region=REGION,
                 cmap=cpt_path, nan_transparent=True)
    gridline_nc.unlink(missing_ok=True)

    # Reconstructed continental polygon outlines (drawn first, in red).
    # Local rotation files (BOTH Clennett + Global 250-0 Ma merged into a
    # single pygplates RotationModel) + local ContinentalPolygons gpml.
    # Silently no-ops if pygplates isn't installed or the input files
    # have moved.
    paleo_continents.plot_continents_on(fig, t,
                                        projection=PROJ, region=REGION,
                                        pen=CONT_POLY_PEN)

    # Reconstructed paleo-coastlines (drawn AFTER the continental polygons,
    # in green, so they sit on top of the red continental-crust outlines).
    # Same local plate model and same merged RotationModel as above -
    # NOT fetched via plate-model-manager.
    paleo_continents.plot_coastlines_on(fig, t,
                                        projection=PROJ, region=REGION,
                                        pen=COASTLINE_PEN)

    # Paleobathymetry is non-positive (depth below sea level); only the low
    # end can overflow the CPT, so add a back (left-side) end-arrow only.
    fig.colorbar(projection=PROJ, region=REGION, cmap=cpt_path,
                 frame=["x+lPaleobathymetry (m)"],
                 position="JBC+w12c/0.3c+o0/1c+h+eb")

    # Time stamp in the top-left, matching the carbonate video offset
    # (-0.15c/-0.1c) so the two videos line up visually if played alongside.
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
    p.add_argument("--cmap", default=DEFAULT_CMAP,
                   help=f"colour palette (default: {DEFAULT_CMAP})")
    p.add_argument("--vmin", type=float, default=DEFAULT_VMIN,
                   help=f"CPT minimum in metres (default: {DEFAULT_VMIN:.0f})")
    p.add_argument("--vmax", type=float, default=DEFAULT_VMAX,
                   help=f"CPT maximum in metres (default: {DEFAULT_VMAX:.0f})")
    p.add_argument("--step", type=float, default=DEFAULT_STEP,
                   help=f"CPT step in metres (default: {DEFAULT_STEP:.0f})")
    args = p.parse_args()

    if args.ages:
        ages = sorted(args.ages, reverse=True)
    else:
        ages = list(range(cfg.MAX_TIME_MA, cfg.MIN_TIME_MA - 1, -args.cadence))

    print(f"\nPaleobathymetry CPT range (fixed): "
          f"{args.vmin:.0f} .. {args.vmax:.0f} m (step {args.step:.0f}, cmap {args.cmap})")
    print(f"Frames: {len(ages)}  (range {ages[0]} -> {ages[-1]} Ma)")

    sub = VIDEO_FRAME_DIR / MODE
    sub.mkdir(exist_ok=True)
    existing = sorted(sub.glob("frame_*.png"))
    if args.force and existing:
        print(f"[{MODE}] --force: wiping {len(existing)} cached frame(s)")
        for f in existing:
            f.unlink()
        existing = []
    elif len(existing) > len(ages):
        print(f"[{MODE}] {len(existing)} stale frames found (need {len(ages)}) - wiping")
        for f in existing:
            f.unlink()
        existing = []
    if existing and len(existing) >= len(ages):
        print(f"[{MODE}] reusing {len(existing)} cached frames "
              f"(pass --force to re-render)")

    for n, t in enumerate(ages):
        frame = sub / f"frame_{n:04d}.png"
        if frame.exists():
            continue
        print(f"[{MODE}] frame {n+1}/{len(ages)} at {t} Ma -> {frame.name}")
        render_frame(t, frame, args.cmap, args.vmin, args.vmax, args.step)

    if args.skip_stitch:
        print(f"[{MODE}] --skip-stitch: frames in {sub}")
        return

    out_video = cfg.OUTPUT_DIR / f"paleobathymetry_{ages[0]}-{ages[-1]}Ma.mp4"
    if out_video.exists():
        out_video.unlink()
    stitch_video(sub / "frame_%04d.png", out_video, args.fps, n_frames=len(ages))
    print(f"  wrote {out_video}")

    if not args.keep_frames:
        shutil.rmtree(sub)


if __name__ == "__main__":
    main()
