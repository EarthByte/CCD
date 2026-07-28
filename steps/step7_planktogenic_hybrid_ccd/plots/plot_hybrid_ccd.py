#!/usr/bin/env python3
"""
Figure — Hybrid CCD (0–170 Ma) with uncertainty envelope + reference curves.

Faithful pyGMT port of ``6_planktogenic_comp/2_ccd_plot.sh`` (same region, pens,
colours and legend). Requires GMT + pyGMT.

Inputs:
    outputs/step7_hybrid/hybrid_ccd_obs_pred_combined.txt   (Age, CCD, ...)
    outputs/step7_hybrid/hybrid_ccd_error_envelope.txt      (closed polygon)
    data/reference_ccd/Boss_Wilkinson_1991_global_CCD_mean.txt
    data/reference_ccd/Global_CCD_Delaney_Boyle.txt

Output:
    outputs/step7_hybrid/hybrid_ccd_with_uncertainty.pdf
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ccdworkflow import config          # noqa: E402
from ccdworkflow.io import read_table, save_pygmt_figure    # noqa: E402

try:
    import pygmt
except ImportError as exc:  # pragma: no cover
    raise SystemExit("pyGMT required: conda install -c conda-forge pygmt") from exc

TMIN, TMAX, YMIN, YMAX = 0, 170, -5000, -2500


def main() -> None:
    db = read_table(config.REFERENCE_CCD / "Global_CCD_Delaney_Boyle.txt", names=["Age", "CCD"])
    db["CCD"] = -db["CCD"]

    pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica", FONT_LABEL="12p,Helvetica",
                 MAP_FRAME_PEN="1p", MAP_TICK_PEN="1p")
    fig = pygmt.Figure()
    fig.basemap(region=[TMIN, TMAX, YMIN, YMAX], projection="X-15c/8c",
                frame=['xa20f5+lAge (Ma)', 'ya500f100+lCCD (m)', "WSen"])
    fig.plot(data=str(config.HYBRID_CCD_ENVELOPE), fill="lightgray", pen="0.25p,lightgray")
    fig.plot(data=str(config.HYBRID_CCD), pen="2.5p,indianred")
    fig.plot(data=str(config.REFERENCE_CCD / "Boss_Wilkinson_1991_global_CCD_mean.txt"),
             pen="2.0p,darkgreen")
    fig.plot(x=db["Age"], y=db["CCD"], pen="2.0p,royalblue3")

    spec = ("S 0.4c - 0.75c - 2.5p,indianred 1c Hybrid CCD\n"
            "S 0.4c - 0.75c - 6p,lightgray 1c Uncertainty\n")
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(spec); spec_path = fh.name
    fig.legend(spec=spec_path, position="jTR")

    png, pdf = save_pygmt_figure(fig, "hybrid_ccd_with_uncertainty", dpi=300, step="step7")
    print(f"[plot] wrote {png} and {pdf}")


if __name__ == "__main__":
    main()
