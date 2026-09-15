#!/usr/bin/env python3
"""Figure 1 (pyGMT): the three regional CCD reconstructions, their inter-basin
range and the area-weighted global mean, 0-52 Ma.

Single-column width for Geology (5.9 cm), so every type size on the page is the
size the reader sees: 7 pt annotations, 8 pt axis labels, 7 pt legend.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt

import paper_gmt as S

CW, OUT = S.CW, S.FIGDIR
W, H = S.W_1COL, 3.6            # cm, frame
AGEMAX = 52

REG = {"Atlantic": ("#1f77b4", CW / "steps/step1_timescale_conversion/outputs/ATL_CCD_GTS2020_1my.txt"),
       "Pacific":  ("#d62728", CW / "steps/step1_timescale_conversion/outputs/PAC_CCD_GTS2020_1my.txt"),
       "Indian":   ("#2ca02c", CW / "steps/step1_timescale_conversion/outputs/IND_CCD_GTS2020_1my.txt")}
GLOBAL = CW / "steps/step3_global_ccd_synthesis/outputs/global_ccd_with_basin_dispersion_0-52Ma.txt"


def rd(p, names):
    return pd.read_csv(p, sep=r"\s+", comment="#", header=None, names=names).dropna()


reg = {k: rd(p, ["age", "ccd"]) for k, (_, p) in REG.items()}
glo = rd(GLOBAL, ["age", "ccd", "lo", "hi"])
for d in list(reg.values()) + [glo]:
    d.drop(d.index[d["age"] > AGEMAX], inplace=True)

# ---- region ----------------------------------------------------------------
lo = min(glo["lo"].min(), min(d["ccd"].min() for d in reg.values()))
hi = max(glo["hi"].max(), max(d["ccd"].max() for d in reg.values()))
Y0, Y1 = np.floor(lo / 250) * 250, np.ceil(hi / 250) * 250
panel = S.Panel(region=(0, AGEMAX, Y0, Y1), width=W, height=H, x_reversed=True)

fig = pygmt.Figure()
pygmt.config(**S.defaults())
S.begin(fig)
fig.basemap(region=list(panel.region), projection=panel.projection,
            frame=["WSne", "xa10f5+lAge (Ma)", "ya500f250+lCCD (m)"])

# inter-basin range first, so the curves sit on top of it
band_x = np.concatenate([glo["age"].to_numpy(), glo["age"].to_numpy()[::-1]])
band_y = np.concatenate([glo["hi"].to_numpy(), glo["lo"].to_numpy()[::-1]])
fig.plot(x=band_x, y=band_y, fill="gray85", pen=None, close=True)

for name, (colour, _) in REG.items():
    d = reg[name]
    fig.plot(x=d["age"], y=d["ccd"], pen=f"0.6p,{colour}")
fig.plot(x=glo["age"], y=glo["ccd"], pen="1.4p,black")

# ---- legend ----------------------------------------------------------------
# Above the frame, in two rows: the panel is 5.9 cm wide and every corner of it
# carries data, so a legend inside the frame sits on the curves it describes.
ENTRIES = [("Atlantic", "#1f77b4", "line", "0.6p"),
           ("Pacific", "#d62728", "line", "0.6p"),
           ("Indian", "#2ca02c", "line", "0.6p"),
           ("Inter-basin range", "gray85", "patch", "0.2p,gray60"),
           ("Global mean", "black", "line", "1.4p")]
legend_box = S.draw_legend(fig, panel, ENTRIES, x0=0.0, y0=H + 0.95, ncol=2)

# ---- checks ----------------------------------------------------------------
bad = []
for name, d in list(reg.items()) + [("Global mean", glo)]:
    if panel.data_in_box(legend_box, d["age"], d["ccd"]):
        bad.append(f"legend / {name} curve")
if panel.data_in_box(legend_box, glo["age"], glo["hi"]):
    bad.append("legend / inter-basin range")
S.report(bad, "legend")
small = S.too_small(S.PT_ANNOT, S.PT_LABEL, S.PT_LEG)
print("  type below the 7 pt journal floor: " + (str(small) if small else "none"))
print(f"  frame {W:.2f} x {H:.2f} cm, legend spans "
      f"{legend_box.x1 - legend_box.x0:.2f} cm of the {W:.2f} cm width")
if legend_box.x1 > W + 0.02:
    print("  WARNING: the legend runs past the panel's right edge")

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Fig1_regional_global_ccd_gmt.{ext}", dpi=600, crop=True)
print("wrote Fig1_regional_global_ccd_gmt")
