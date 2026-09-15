#!/usr/bin/env python3
"""Figure 4 (pyGMT), two-column width for Geology (12.28 cm).

(A) the global CCD against the solid-Earth carbon fluxes recomputed from it - the
    total the attribution uses, its four components, and the one sink the model
    computes, drawn with the opposite sign;
(B) zero-lag correlation of sea level and of each flux with the CCD;
(C) variance of the CCD that sea level and total degassing each explain alone.

Reads attribution_matched_series.csv and component_correlations.csv, both written
by attribution_closure_analysis.py, so the panels and the reported numbers cannot
drift apart.
"""
import importlib.util as _ilu

import numpy as np
import pandas as pd
import pygmt

import paper_gmt as S

OUT = S.FIGDIR
W = S.W_2COL
HA = 4.4                     # cm, panel A frame
H2 = 3.8                     # cm, panels B and C
ROW_GAP = 1.30               # cm between panel A's age label and the row below
LEG_ROW = 1.20               # cm above panel A for its legend, clear of the
                             # event labels below it
AGEMAX = 170

D = pd.read_csv(S._HERE / "attribution_matched_series.csv")
_CCFILE = S._HERE / "component_correlations.csv"
if not _CCFILE.exists():
    raise SystemExit(f"{_CCFILE.name} not found - run attribution_closure_analysis.py first")
CC = pd.read_csv(_CCFILE)
age = D["age"].to_numpy()
ccd = -D["depth"].to_numpy()

# One colour per series, wherever that series appears: the curves in (A), the markers
# in (B) and the bars in (C).
COL = {k: None for k in ()} or {"Sea level": "#009e73", "Mid-ocean ridge": "#0072b2", "Rift": "#9467bd",
       "Carbonate platform": "#8c564b", "Arc (subduction)": "#d1495b",
       "Total outflux": "#3d3d3d", "Seafloor weathering": "#e69f00"}

COL = {k: S.rgb(v) for k, v in COL.items()}      # GMT pens want r/g/b, not hex

_espec = _ilu.spec_from_file_location("paper_events", S._HERE / "paper_events.py")
_ev = _ilu.module_from_spec(_espec)
_espec.loader.exec_module(_ev)

fig = pygmt.Figure()
pygmt.config(**S.defaults())
S.begin(fig)
checks = []

# ============================ row 2: (B) and (C) =============================
# Drawn first because they sit at the bottom of the page and GMT builds upward.
BX, BW_ = 2.45, 3.45         # (B): its category labels live in the 2.45 cm to its left
CX, CW_ = 8.30, 3.75         # (C)
pB = S.Panel(region=(-1.0, 1.0, -0.7, len(CC) - 0.3), width=BW_, height=H2, origin=(BX, 0.0))
pC = S.Panel(region=(-0.6, 1.6, 0.0, 1.0), width=CW_, height=H2, origin=(CX, 0.0))

# ---- (B) correlations -------------------------------------------------------
fig.shift_origin(xshift=f"{BX}c")
RB, PB = list(pB.region), pB.projection
fig.basemap(region=RB, projection=PB, frame=["wSrt", "xa0.5f0.25+lCorrelation with CCD, 0-52 Ma", "ya1"])
fig.plot(x=[0, 0], y=[RB[2], RB[3]], pen="0.5p,gray35", region=RB, projection=PB)
_cc = CC.iloc[::-1].reset_index(drop=True)      # first row of the file at the top
cat_boxes = []
for i, row in _cc.iterrows():
    colour = COL.get(row["component"], "gray40")
    # Squares mark the carbon sink, circles the sources and sea level; a sink enters
    # with its sign reversed, as its contribution to CO2, so positive always means the
    # sign a CO2 control requires. Open symbols are computed from the CCD itself.
    style = "s0.20c" if bool(row.get("is_sink", False)) else "c0.22c"
    fig.plot(x=[row["r"]], y=[i], style=style, pen=f"0.9p,{colour}",
             fill=colour if row["ccd_independent"] else "white", region=RB, projection=PB)
    fig.text(x=RB[0], y=i, text=row["component"], justify="RM", offset="-0.12c/0c",
             font=f"{S.PT_ANNOT}p,{S.FONT},black", no_clip=True, region=RB, projection=PB)
    cat_boxes.append(pB.label_box(row["component"], row["component"], RB[0], i,
                                  S.PT_ANNOT, justify="RM", dx=-0.12))
checks += S.collisions(cat_boxes)
if min(b.x0 for b in cat_boxes) < -BX + 0.05:
    checks.append("category labels run off the left edge of the figure")
fig.text(x=-0.62, y=H2 + 0.28, text="B", justify="LB", font=f"{S.PT_TAG}p,{S.FONT}-Bold,black",
         no_clip=True, **S.cm_frame(pB, top=0.9))
fig.shift_origin(xshift=f"-{BX}c")

# ---- (C) variance explained -------------------------------------------------
def z(a):
    a = np.asarray(a, float)
    return (a - a.mean()) / a.std()


sub = D[D["age"] <= 52]
Y = z(-sub["depth"])            # CCD in the sense the figures plot it


def r2_of(*cols):
    X = np.column_stack([np.ones(len(Y))] + list(cols))
    b, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    return 1 - np.sum((Y - X @ b) ** 2) / np.sum((Y - Y.mean()) ** 2)


r2_sl = r2_of(z(sub["sl"]))
r2_dg = r2_of(z(sub["gross_outflux"]))

fig.shift_origin(xshift=f"{CX}c")
RC, PC = list(pC.region), pC.projection
# Acidification by added CO2 requires a POSITIVE coefficient on this axis; degassing
# takes a negative one, so its variance is explained in the wrong direction and cannot
# support a CO2 control however large it is. The hatch and the label say so.
CLABEL = "Variance explained (R@+2@+)"
# Lower-case s: the bottom axis is drawn but not annotated, because the two bars are
# named underneath it in words.
fig.basemap(region=RC, projection=PC, frame=["Wsrt", f"ya0.2f0.1+l{CLABEL}"])
for x, val, name, hatched in ((0, r2_sl, "Sea level", False), (1, r2_dg, "Total outflux", True)):
    # White hatch on the series colour, so the bar still reads as total outflux.
    fill = COL[name] if not hatched else f"p8+fwhite+b{COL[name]}"
    fig.plot(x=[x], y=[val], style="b0.55u+b0", fill=fill, pen="0.5p,black",
             region=RC, projection=PC)
    fig.text(x=x, y=val + 0.02, text=f"{val:.2f}", justify="CB",
             font=f"{S.PT_ANNOT}p,{S.FONT},gray15", region=RC, projection=PC)
    fig.text(x=x, y=0.0, text=name, justify="CT", offset="0c/-0.12c",
             font=f"{S.PT_ANNOT}p,{S.FONT},black", no_clip=True, region=RC, projection=PC)
fig.text(x=1, y=r2_dg + 0.10, text="wrong sign", justify="CB",
         font=f"{S.PT_ANNOT}p,{S.FONT},gray15", region=RC, projection=PC)
_lab_len = S.text_width_cm(CLABEL.replace("@+", ""), S.PT_LABEL)
if _lab_len > H2:
    checks.append(f"(C) y-label is {_lab_len:.1f} cm long on a {H2:.1f} cm panel")
fig.text(x=-0.55, y=H2 + 0.28, text="C", justify="LB", font=f"{S.PT_TAG}p,{S.FONT}-Bold,black",
         no_clip=True, **S.cm_frame(pC, top=0.9))
fig.shift_origin(xshift=f"-{CX}c")

# ============================ panel A ========================================
fig.shift_origin(yshift=f"{H2 + ROW_GAP}c")
pA = S.Panel(region=(0, AGEMAX, -5100, -2500), width=W, height=HA, x_reversed=True)
pF = S.Panel(region=(0, AGEMAX, -60, 175), width=W, height=HA, x_reversed=True)
RA, PA = list(pA.region), pA.projection
fig.basemap(region=RA, projection=PA,
            frame=["WSn", "xa20f10+lAge (Ma)", "ya500f250+lCCD (m)"])
# 2.6 mm off the frame rather than 1.6: at this panel width the labels sat on the
# topmost depth annotation.
_event_boxes = _ev.draw_events_gmt(fig, pA, pt=S.PT_LEG, dy=0.26)
fig.plot(x=age, y=ccd, pen="1.4p,black", region=RA, projection=PA)

RF = list(pF.region)
fig.basemap(region=RF, projection=PA, frame=["E", "ya50f25+lCarbon flux (Mt C yr@+-1@+)"])
fig.plot(x=[0, AGEMAX], y=[0, 0], pen="0.4p,gray55", region=RF, projection=PA)
SERIES = [("gross_outflux", "Total outflux", "1.2p,%s,4_1.5:0", 1),
          ("arc_subduction", "Arc (subduction)", "0.8p,%s", 1),
          ("MOR_ridge", "Mid-ocean ridge", "0.8p,%s", 1),
          ("rift", "Rift", "0.6p,%s", 1),
          ("carb_platform", "Carbonate platform", "0.6p,%s", 1),
          ("crust_sink", "Seafloor weathering", "0.7p,%s,3_1.2_0.6_1.2:0", -1)]
for col, name, pen, sign in SERIES:
    fig.plot(x=age, y=sign * D[col].to_numpy(), pen=pen % COL[name], region=RF, projection=PA)

ENTRIES_A = [("Global CCD (this study)", "black", "line", "1.4p"),
             ("Total outflux (model predictor)", COL["Total outflux"], "line", "1.2p,4_1.5:0"),
             ("Arc outgassing (pelagic carbonates)", COL["Arc (subduction)"], "line", "0.8p"),
             ("Mid-ocean ridge outgassing", COL["Mid-ocean ridge"], "line", "0.8p"),
             ("Rift outgassing", COL["Rift"], "line", "0.6p"),
             ("Arc outgassing (carbonate platforms)", COL["Carbonate platform"], "line", "0.6p"),
             ("Seafloor weathering (carbon uptake)", COL["Seafloor weathering"],
              "line", "0.7p,3_1.2_0.6_1.2:0")]
# Above the panel: seven entries over a curve that fills the frame had nowhere to sit
# inside it. The row clears the event labels, which sit just above the frame.
legA = S.draw_legend(fig, pA, ENTRIES_A, x0=0.10, y0=HA + LEG_ROW + 0.80, ncol=2, col_gap=0.30)
if legA.x1 > W + 0.05:
    checks.append(f"legend A is {legA.x1 - legA.x0:.1f} cm wide on a {W:.1f} cm figure")
# The legend's bottom row and the event labels share the strip above the panel, so the
# clearance between them is measured rather than judged by eye.
_gap = legA.y0 - max(b.y1 for b in _event_boxes)
print(f"  legend A clears the event labels by {_gap * 10:.1f} mm")
if _gap < 0.20:
    checks.append(f"legend A sits {_gap * 10:.1f} mm above the event labels")

fig.text(x=-0.62, y=HA + 0.28, text="A", justify="LB", font=f"{S.PT_TAG}p,{S.FONT}-Bold,black",
         no_clip=True, **S.cm_frame(pA, top=0.9))

# ============================ checks and output ==============================
_pmin = float(CC["p_ar1"].min())
S.report(checks, "layout")
small = S.too_small(S.PT_ANNOT, S.PT_LABEL, S.PT_LEG)
print("  type below the 7 pt journal floor: " + (str(small) if small else "none"))
print(f"  caption facts: positive = flux rises as the CCD shoals; open symbols are "
      f"computed from the CCD; all p > {np.floor(_pmin * 10) / 10:.1f}")
print(f"  R2 alone: sea level {r2_sl:.2f}, total degassing {r2_dg:.2f} (wrong sign)")

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Fig4_combined_gmt.{ext}", dpi=600, crop=True)
print("wrote Fig4_combined_gmt")
