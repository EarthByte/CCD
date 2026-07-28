#!/usr/bin/env python3
"""
=============================================================================
06_delta_panel.py  -  2 x 2 global delta maps at the four picked times
=============================================================================

Reads the four picked times from output/stats/picked_times.txt (written by
02_difference_stats.py) and lays out the DM2026 - BW1991 carbonate-
thickness difference field at each in a 2 x 2 grid, all sharing the same
polar cpt and one common colorbar across the bottom.

Same global Winkel-Tripel projection, same Alfonso 2024 reconstructed
coastlines, and the same +e end-arrow colorbar style as the delta video,
so the panel reads as four still frames from the same animation.

Output:
    output/figures/06_delta_panel.png
    output/figures/06_delta_panel.pdf

Run AFTER 02_difference_stats.py.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pygmt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as cfg
import paleo_coastlines


PICKED_TIMES_PATH = cfg.OUTPUT_DIR / "stats" / "picked_times.txt"

# Per-panel projection width (so the 2x2 fits a single-column figure).
PROJ_WIDTH_CM = 9
PROJ = f"R0/{PROJ_WIDTH_CM}c"
# REGION = "d" (-180/180/-90/90) matches the native longitude convention of
# the carbonate grids and avoids the "Longitude range too small" warning
# that grdimage emits when the plot region's lon origin differs from the
# grid's. See cfg.to_gridline_nc() for the matching registration normaliser.
REGION = "d"

CMAP = "polar"
DELTA_CAP_PCT = 99


def load_picked_times() -> list[int]:
    if not PICKED_TIMES_PATH.exists():
        raise SystemExit(f"{PICKED_TIMES_PATH} missing -- run "
                         "02_difference_stats.py first.")
    picks = []
    with open(PICKED_TIMES_PATH) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("time_Ma"):
                continue
            picks.append(int(line.split(",")[0]))
    if len(picks) != 4:
        print(f"  [warn] expected 4 picked times, got {len(picks)}: {picks}")
    return picks


def compute_delta_range(times) -> tuple[float, float, float]:
    """Symmetric range capped at the p99 of |delta| across the four panels."""
    abs_max = 0.0
    for t in times:
        d = cfg.diff_grid(t).values
        d = d[np.isfinite(d)]
        if d.size == 0:
            continue
        c = float(np.percentile(np.abs(d), DELTA_CAP_PCT))
        if c > abs_max:
            abs_max = c
    for step in (5, 10, 20, 25, 50, 100, 200):
        if abs_max / step <= 25:
            break
    cap = float(np.ceil(abs_max / step) * step)
    return -cap, cap, float(step)


def main() -> int:
    picks = load_picked_times()[:4]
    print(f"  picked times for 2x2 panel: {picks}")

    lo, hi, step = compute_delta_range(picks)
    print(f"  shared delta range across panels: +/- {hi:.0f} m (step {step:.0f})")

    cpt_path = "_delta_panel.cpt"
    pygmt.makecpt(cmap=CMAP, series=[lo, hi, step], continuous=True,
                  background=False, output=cpt_path)

    fig = pygmt.Figure()
    pygmt.config(MAP_FRAME_TYPE="plain",
                 FONT_ANNOT_PRIMARY="8p,Helvetica,black",
                 FONT_LABEL="10p,Helvetica,black",
                 FONT_TITLE="11p,Helvetica-Bold,black",
                 COLOR_NAN="220/220/220",
                 MAP_TITLE_OFFSET="-2p")

    with fig.subplot(
        nrows=2, ncols=2,
        figsize=(PROJ_WIDTH_CM * 2 + 0.5, PROJ_WIDTH_CM * 2 * 0.55),
        margins=["0.3c", "0.3c"],
        sharex="b", sharey="l",
    ):
        for i, t in enumerate(picks):
            with fig.set_panel(panel=i):
                da = cfg.diff_grid(t)
                fig.basemap(region=REGION, projection=PROJ,
                            frame=["af", f'+t"{t} Ma"'])
                # Force gridline registration on the way into grdimage to
                # silence the "Longitude range too small" warning regardless
                # of how the upstream grid was produced.
                gridline_nc = cfg.to_gridline_nc(da)
                fig.grdimage(grid=str(gridline_nc),
                             projection=PROJ, region=REGION,
                             cmap=cpt_path, nan_transparent=True)
                gridline_nc.unlink(missing_ok=True)
                paleo_coastlines.plot_on(fig, t,
                                         projection=PROJ, region=REGION,
                                         pen="0.35p,gray30")
                # Panel tag in the top-left.
                tag = "(a)(b)(c)(d)"[i*3:i*3 + 3]
                fig.text(text=tag, font="11p,Helvetica-Bold,black",
                         position="TL", justify="TL",
                         offset="-0.1c/0.6c", no_clip=True,
                         region=REGION, projection=PROJ)

    # Single shared colorbar across the bottom.
    fig.colorbar(cmap=cpt_path,
                 frame=['x+l"thickness anomaly (m), DM2026 - BW1991"'],
                 position="JBC+w14c/0.3c+o0/1.5c+h+e")

    out_png = cfg.OUTPUT_DIR / "figures" / "06_delta_panel.png"
    out_pdf = cfg.OUTPUT_DIR / "figures" / "06_delta_panel.pdf"
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    print(f"  wrote {out_png}")
    print(f"  wrote {out_pdf}")

    if Path(cpt_path).exists():
        Path(cpt_path).unlink()
    return 0


if __name__ == "__main__":
    sys.exit(main())
