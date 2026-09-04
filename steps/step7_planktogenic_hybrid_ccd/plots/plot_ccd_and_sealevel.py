#!/usr/bin/env python3
"""
Figure — CCD (left axis) + long-term sea level (right axis), 0–170 Ma.

Faithful pyGMT port of ``6_planktogenic_comp/3_ccd_and_sl_plot.sh``: primary CCD
basemap on the left, an overlaid secondary basemap for the sea-level axis on the
right, identical pens/colours and legend. Requires GMT + pyGMT.

Inputs:
    outputs/step6_regression/predicted_ccd_sl_0-205Ma.txt       (sea-level CCD)
    outputs/step7_hybrid/hybrid_ccd_obs_pred_combined.txt       (hybrid CCD)
    outputs/step7_hybrid/hybrid_ccd_error_envelope.txt          (envelope)
    outputs/step4_sealevel_envelope/sea_level_quantile_envelope_0-205Ma.txt (this study SL)
    data/reference_ccd/Boss_Wilkinson_1991_global_CCD_mean.txt
    data/reference_ccd/Global_CCD_Delaney_Boyle.txt
    steps/step1_timescale_conversion/outputs/Haq87_longterm_GTS2020.txt

Output:
    outputs/step7_hybrid/CCD_SeaLevel.pdf  (+ .png)
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ccdworkflow import config          # noqa: E402
from ccdworkflow.io import read_table, read_xy, save_pygmt_figure   # noqa: E402

try:
    import pygmt
except ImportError as exc:  # pragma: no cover
    raise SystemExit("pyGMT required: conda install -c conda-forge pygmt") from exc

TMIN, TMAX = 0, 170
CCD_MIN, CCD_MAX = -5300, -2500
SL_MIN, SL_MAX = 0, 300


def _map_sl(df):
    scale = (CCD_MAX - CCD_MIN) / (SL_MAX - SL_MIN)
    offset = CCD_MIN - SL_MIN * scale
    return df.iloc[:, 0].to_numpy(), df.iloc[:, 1].to_numpy() * scale + offset


def main() -> None:
    db = read_table(config.REFERENCE_CCD / "Global_CCD_Delaney_Boyle.txt", names=["Age", "CCD"])
    db["CCD"] = -db["CCD"]
    sl = read_xy(config.SL_ENVELOPE_FULL, names=["Age", "SL"])
    haq = read_xy(config.GTS2020_HAQ_LONGTERM, names=["Age", "SL"])
    sl_x, sl_y = _map_sl(sl)
    haq_x, haq_y = _map_sl(haq)

    pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica", FONT_LABEL="12p,Helvetica",
                 MAP_FRAME_PEN="1p", MAP_TICK_PEN="1p")
    fig = pygmt.Figure()
    fig.basemap(region=[TMIN, TMAX, CCD_MIN, CCD_MAX], projection="X-15c/8c",
                frame=['xa20f5+lAge (Ma)', 'ya500f100+lCCD (m)', "WSn"])
    fig.plot(data=str(config.HYBRID_CCD_ENVELOPE), fill="lightgray", pen="0.25p,lightgray")
    fig.plot(data=str(config.PREDICTED_CCD_SL), pen="2.0p,cadetblue")
    fig.plot(data=str(config.HYBRID_CCD), pen="3.0p,mediumturquoise")
    fig.plot(data=str(config.REFERENCE_CCD / "Boss_Wilkinson_1991_global_CCD_mean.txt"),
             pen="2.0p,darkorange2")
    fig.plot(x=db["Age"], y=db["CCD"], pen="2.0p,magenta4")
    fig.plot(x=sl_x, y=sl_y, pen="1.25p,royalblue4")
    fig.plot(x=haq_x, y=haq_y, pen="1.25p,skyblue3")

    # secondary right axis in sea-level units
    fig.basemap(region=[TMIN, TMAX, SL_MIN, SL_MAX], projection="X-15c/8c",
                frame=['ya50f10+lSea level (m)', "E"])

    spec = (
        "H 12p CCD\n"
        "S 0.4c - 0.75c - 2.0p,cadetblue 1c Sealevel CCD\n"
        "S 0.4c - 0.75c - 2.5p,mediumturquoise 1c Hybrid CCD\n"
        "S 0.4c - 0.75c - 6p,lightgray 1c CCD uncertainty\n"
        "S 0.4c - 0.75c - 2.5p,darkorange2 1c Boss & Wilkinson (1991)\n"
        "S 0.4c - 0.75c - 2.5p,magenta4 1c Delaney & Boyle\n"
        "G 0.15c\n"
        "H 12p Sea Level\n"
        "S 0.4c - 0.75c - 1.25p,royalblue4 1c This study\n"
        "S 0.4c - 0.75c - 1.25p,skyblue3 1c Haq (1987)\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(spec); spec_path = fh.name
    fig.legend(spec=spec_path, position="jBC+o1.4c/1.5c+w8c/2.6c")

    png, pdf = save_pygmt_figure(fig, "CCD_SeaLevel", dpi=300, step="step7")
    print(f"[plot] wrote {png} and {pdf}")


if __name__ == "__main__":
    main()
