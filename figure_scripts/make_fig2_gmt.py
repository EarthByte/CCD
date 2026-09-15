#!/usr/bin/env python3
"""Figure 2 (pyGMT), two-column width for Geology (12.28 cm).

(A) the long-term sea-level record with the post-52 Ma calibration of the CCD
    against sea level as an inset;
(B) the global CCD since 170 Ma with its uncertainty, the earlier compilations
    and the sea-level records on a right-hand axis.
"""
import importlib.util as _ilu

import numpy as np
import pandas as pd
import pygmt

import paper_gmt as S

CW, OUT = S.CW, S.FIGDIR
W = S.W_2COL
HA = HB = 4.6          # cm, frame height of each panel
GAP = 1.55             # cm between the frames: A's age label and B's event labels
AGEMAX = 170
MHT = 66.6             # Miller/Haq transition

COL = dict(short="gray72", haq="#2b5fa6", env="#c0392b",
           pred="#5f9ea0", hyb="#48d1cc", band="#d3d3d3",
           bw="#ee7600", db="#8b008b")


def rd(p, names):
    return pd.read_csv(p, sep=r"\s+", engine="python", comment="#", header=None, names=names)


sl = rd(CW / "steps/step1_timescale_conversion/outputs/sealevel_shortterm_hybrid_GTS2020.txt",
        ["age", "sl"]).sort_values("age")
haq = rd(CW / "data/sealevel/Haq87_Longterm_v3.txt", ["age", "sl"]).sort_values("age")
env = rd(CW / "steps/step4_sealevel_envelope/outputs/sea_level_quantile_envelope_0-205Ma.txt",
         ["age", "sl"]).sort_values("age")
mt = pd.read_csv(CW / "steps/step6_sealevel_ccd_regression/outputs/diagnostics/matched_timeseries.csv")
rr = pd.read_csv(CW / "steps/step6_sealevel_ccd_regression/outputs/diagnostics/regression_results.csv")
rma = rr[rr["method"] == "RMA"].iloc[0]
_bandf = CW / "steps/step6_sealevel_ccd_regression/outputs/diagnostics/bootstrap_band_rma.csv"
band = pd.read_csv(_bandf) if _bandf.exists() else None
hyb = rd(OUT / "CCD_hybrid_DM2026.txt", ["age", "ccd", "lo", "hi"]).dropna()
pred = rd(CW / "steps/step6_sealevel_ccd_regression/outputs/predicted_ccd_sl_0-205Ma.txt",
          ["age", "ccd", "lo", "hi"]).dropna()
bw = rd(CW / "data/reference_ccd/Boss_Wilkinson_1991_global_CCD_mean.txt",
        ["age", "ccd", "mn", "mx", "a2", "mx2", "a3", "mn3"]).dropna(subset=["age", "ccd"])
db = rd(CW / "data/reference_ccd/Global_CCD_Delaney_Boyle.txt", ["age", "ccd"]).dropna()
# This file is depth positive downward (3298-4491 m) while the panel plots the CCD
# negative downward, so without the flip the curve lands off the top of the frame -
# which is where it was in the matplotlib figure: legended but never visible.
db["ccd"] = -db["ccd"]

_espec = _ilu.spec_from_file_location("paper_events", S._HERE / "paper_events.py")
_ev = _ilu.module_from_spec(_espec)
_espec.loader.exec_module(_ev)

fig = pygmt.Figure()
pygmt.config(**S.defaults())
S.begin(fig)
checks = []
LEG_ROW = 0.80      # cm, the strip under panel B that holds its legend
AXIS_ROW = 1.00     # cm, panel B's age annotations and label
fig.shift_origin(yshift=f"{LEG_ROW + AXIS_ROW}c")

# ============================ panel B (drawn first, at the origin) ===========
pb = S.Panel(region=(0, AGEMAX, -5750, -2500), width=W, height=HB, x_reversed=True)
pb_sl = S.Panel(region=(0, AGEMAX, 0, 300), width=W, height=HB, x_reversed=True)
RB, PB = list(pb.region), pb.projection
fig.basemap(region=RB, projection=PB,
            frame=["WSn", "xa20f10+lAge (Ma)", "ya500f250+lCCD (m)"])

fig.plot(x=np.concatenate([hyb["age"], hyb["age"][::-1]]),
         y=np.concatenate([hyb["hi"], hyb["lo"][::-1]]),
         fill=COL["band"], close=True, region=RB, projection=PB)
_ev.draw_events_gmt(fig, pb, pt=S.PT_LEG)
for d, c, w in ((pred, COL["pred"], "0.7p"), (bw, COL["bw"], "0.7p"),
                (db, COL["db"], "0.7p"), (hyb, COL["hyb"], "1.3p")):
    fig.plot(x=d["age"], y=d["ccd"], pen=f"{w},{c}", region=RB, projection=PB)

# right-hand sea-level axis: same box, its own region
RS = list(pb_sl.region)
fig.basemap(region=RS, projection=PB, frame=["E", "ya100f50+lSea level (m)"])
fig.plot(x=env["age"], y=env["sl"], pen=f"0.6p,{COL['env']}", region=RS, projection=PB)
fig.plot(x=haq["age"], y=haq["sl"], pen=f"0.6p,{COL['haq']}", region=RS, projection=PB)

ENTRIES_B = [("Sea-level CCD", COL["pred"], "line", "0.7p"),
             ("Hybrid CCD", COL["hyb"], "line", "1.3p"),
             ("CCD uncertainty", COL["band"], "patch", "0.2p,gray60"),
             ("Boss & Wilkinson (1991)", COL["bw"], "line", "0.7p"),
             ("Delaney & Boyle (1988)", COL["db"], "line", "0.7p"),
             ("Sea level: this study", COL["env"], "line", "0.6p"),
             ("Sea level: Haq (1987)", COL["haq"], "line", "0.6p")]
# Below the panel, under the age label: the sea-level curve reaches the bottom right
# corner of the frame and the CCD curves the bottom left, so no strip inside the frame
# is free for seven entries. Four columns keep it two rows deep.
fig.shift_origin(yshift=f"-{LEG_ROW + AXIS_ROW}c")
_legpanel = S.Panel(region=(0, 1, 0, 1), width=W, height=LEG_ROW)
legB = S.draw_legend(fig, _legpanel, ENTRIES_B, x0=0.15, y0=LEG_ROW - 0.04,
                     ncol=4, col_gap=0.18)
fig.shift_origin(yshift=f"{LEG_ROW + AXIS_ROW}c")

fig.text(x=-0.55, y=HB + 0.10, text="B", justify="LB", font=f"{S.PT_TAG}p,{S.FONT}-Bold,black",
         no_clip=True, **S.cm_frame(pb, top=0.9))

# ============================ panel A ========================================
fig.shift_origin(yshift=f"{HB + GAP}c")
pa = S.Panel(region=(0, AGEMAX, -60, 275), width=W, height=HA, x_reversed=True)
RA, PA = list(pa.region), pa.projection
fig.basemap(region=RA, projection=PA,
            frame=["WSne", "xa20f10+lAge (Ma)", "ya50f25+lSea level (m)"])
fig.plot(x=[0, AGEMAX], y=[0, 0], pen="0.5p,black,-", region=RA, projection=PA)
fig.plot(x=sl["age"], y=sl["sl"], pen=f"0.3p,{COL['short']}", region=RA, projection=PA)
fig.plot(x=haq["age"], y=haq["sl"], pen=f"1.0p,{COL['haq']}", region=RA, projection=PA)
fig.plot(x=env["age"], y=env["sl"], pen=f"1.2p,{COL['env']}", region=RA, projection=PA)
fig.plot(x=[MHT, MHT], y=[-60, 275], pen="0.6p,gray35,-", region=RA, projection=PA)

# Which record the short-term curve, and hence the red envelope, comes from on
# each side of the transition.
SRC_PT = 8.0
srcs = [("HaqHybrid", MHT + 3.0, "RT"), ("Miller2024", MHT - 3.0, "LT")]
for text, age, just in srcs:
    fig.text(x=age, y=268, text=text, justify=just, font=f"{SRC_PT}p,{S.FONT}-Bold,{COL['env']}",
             region=RA, projection=PA)
src_boxes = [pa.label_box(t, t, a, 268, SRC_PT, justify=j, bold=True) for t, a, j in srcs]

ENTRIES_A = [("Short-term record", COL["short"], "line", "0.3p"),
             ("Haq (1987)", COL["haq"], "line", "1.0p"),
             ("Peak-following envelope", COL["env"], "line", "1.2p")]
legA = S.draw_legend(fig, pa, ENTRIES_A, x0=0.20, y0=HA - 0.18, ncol=1)
for name, d in (("short-term record", sl), ("Haq", haq), ("envelope", env)):
    if pa.data_in_box(legA, d["age"], d["sl"]):
        checks.append(f"legend A / {name}")
checks += S.collisions([legA] + src_boxes)

# ---- inset: the RMA calibration --------------------------------------------
IX, IY, IW, IH = 4.40, 0.75, 4.20, 1.95          # cm within panel A
ins = S.Panel(region=(-25, 150, 3250, 4750), width=IW, height=IH,
              origin=(IX, IY), y_reversed=True)
inset_box = S.Box("inset", IX - 0.95, IY - 0.75, IX + IW + 0.10, IY + IH + 0.10)
for name, d in (("short-term record", sl), ("Haq", haq), ("envelope", env)):
    if pa.data_in_box(inset_box, d["age"], d["sl"]):
        checks.append(f"inset / {name}")
checks += S.collisions([inset_box, legA] + src_boxes)

# The backing covers the inset's own annotations and axis labels as well as its frame,
# so the sea-level curves do not run through the numbers.
_cfa = S.cm_frame(pa)
fig.plot(x=[inset_box.x0, inset_box.x1, inset_box.x1, inset_box.x0],
         y=[inset_box.y0, inset_box.y0, inset_box.y1, inset_box.y1],
         fill="white", pen="0.3p,gray60", close=True, **_cfa)

fig.shift_origin(xshift=f"{IX}c", yshift=f"{IY}c")
RI, PI = list(ins.region), ins.projection
fig.plot(x=[RI[0], RI[1], RI[1], RI[0]], y=[RI[2], RI[2], RI[3], RI[3]],
         fill="white", pen="0.4p,gray40", close=True, region=RI, projection=PI)
if band is not None:
    fig.plot(x=np.concatenate([band["Sea_level_m"], band["Sea_level_m"][::-1]]),
             y=np.concatenate([band["CCD_p2.5"], band["CCD_p97.5"][::-1]]),
             fill="gray70", close=True, region=RI, projection=PI)
fig.plot(x=mt["Sea_level_m"], y=mt["CCD_m_pos_down"], style="c0.09c", pen="0.4p,black",
         region=RI, projection=PI)
_xx = np.linspace(mt["Sea_level_m"].min(), mt["Sea_level_m"].max(), 100)
fig.plot(x=_xx, y=rma["slope"] * _xx + rma["intercept"], pen="1.0p,black",
         region=RI, projection=PI)
fig.basemap(region=RI, projection=PI,
            frame=["WSne", "xa50f25+lSea level (m)", "ya500f250+lCCD (m)"])
# The fit is on CCD positive downward; every correlation in the paper takes the CCD
# in the sense the figures plot it, so no signed r is quoted here - the slope's sign
# says high sea level goes with a shallow CCD and R@+2@+ carries the strength.
EQ_PT = 7.0
eq1 = f"CCD = {rma['intercept']:.0f} - {abs(rma['slope']):.2f} x SL"
eq2 = f"R@+2@+ = {rma['r2']:.2f}"
fig.text(x=RI[0] + 4, y=RI[2] + 60, text=eq1, justify="LT", font=f"{EQ_PT}p,{S.FONT},black",
         fill="white@20", pen="0.3p,gray50", clearance="0.06c/0.04c", region=RI, projection=PI)
fig.text(x=RI[0] + 4, y=RI[2] + 260, text=eq2, justify="LT", font=f"{EQ_PT}p,{S.FONT},black",
         fill="white@20", pen="0.3p,gray50", clearance="0.06c/0.04c", region=RI, projection=PI)
fig.shift_origin(xshift=f"-{IX}c", yshift=f"-{IY}c")

fig.text(x=-0.55, y=HA + 0.10, text="A", justify="LB", font=f"{S.PT_TAG}p,{S.FONT}-Bold,black",
         no_clip=True, **S.cm_frame(pa, top=0.9))

# ============================ checks and output ==============================
S.report(checks, "layout")
small = S.too_small(S.PT_ANNOT, S.PT_LABEL, S.PT_LEG, SRC_PT, EQ_PT)
print("  type below the 7 pt journal floor: " + (str(small) if small else "none"))
print(f"  panels {W:.2f} x {HA:.2f} cm, legend B {legB.x1 - legB.x0:.2f} cm wide")
if legB.x1 > W:
    print("  WARNING: legend B runs past the panel's right edge")

for ext in ("png", "pdf"):
    fig.savefig(OUT / f"Fig2_combined_gmt.{ext}", dpi=600, crop=True)
print("wrote Fig2_combined_gmt")
