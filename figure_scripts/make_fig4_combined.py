#!/usr/bin/env python3
"""Fig. 4 (combined):
   (a) global CCD against the recomputed solid-Earth carbon outflux, by component and
       in total, the total being the predictor the two-predictor model actually uses;
   (b) incremental variance of the two predictors, 0-52 Ma;
   (c) zero-lag correlation of each outflux component with the CCD, 0-52 Ma.
Reads attribution_matched_series.csv and component_correlations.csv, both written by
attribution_closure_analysis.py, so the panels and the reported numbers cannot drift."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- path resolution -------------------------------------------------------
# Works both in the working tree (this script in Paper/, the workflow in a
# sibling CCD_workflow_clean/) and in the public repo (this script in
# figure_scripts/, the workflow at the repo root). CW = workflow root,
# FIGDIR = where figures and the small figure-input curves live.
_HERE = Path(__file__).resolve().parent
def _workflow_root():
    for _p in (_HERE.parent / "CCD_workflow_clean", _HERE.parent, _HERE):
        if (_p / "steps").is_dir() and (_p / "ccdworkflow").is_dir():
            return _p
    raise SystemExit("cannot locate the CCD workflow root (expected steps/ + ccdworkflow/)")
CW = _workflow_root()
FIGDIR = _HERE / "Figures" if (_HERE / "Figures").is_dir() else CW / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)
FIG = str(FIGDIR)
D = pd.read_csv(_HERE/"attribution_matched_series.csv")
_CCFILE = _HERE/"component_correlations.csv"
if not _CCFILE.exists():
    raise SystemExit(f"{_CCFILE.name} not found - run attribution_closure_analysis.py first")
CC = pd.read_csv(_CCFILE)
age=D["age"].to_numpy(); ccd=-D["depth"].to_numpy()

# One colour per outflux component, used by both panel (a) and panel (c).
COL = {"Mid-ocean ridge":"#0072b2", "Rift":"#9467bd", "Carbonate platform":"#8c564b",
       "Arc (subduction)":"#d1495b", "Total outflux":"#3d3d3d"}
PT_TICK, PT_LABEL, PT_LEG = 9.5, 10, 9

# Panel (a) spans the width; (b) and (c) share the row below it. The lag panel is gone:
# carbonate compensation restores ocean saturation in about 5,000-10,000 years, whether
# the perturbation is a sea-level-driven shift in carbonate burial or an added
# solid-Earth CO2 input, which is 100 to 200 times shorter than the 1 Myr sampling. A
# lead or lag of several Myr is not resolvable here, so only zero lag is shown.
fig=plt.figure(figsize=(7.48,8.6))
gs=fig.add_gridspec(2,2,height_ratios=[1.34,0.95],hspace=0.34,wspace=0.62,
                    top=0.835,bottom=0.075,left=0.095,right=0.895)
axa=fig.add_subplot(gs[0,:]); axb1=fig.add_subplot(gs[1,0]); axc=fig.add_subplot(gs[1,1])

# ---- (a) CCD against the carbon outflux -------------------------------------
axa.plot(age,ccd,color="black",lw=2.4,label="Global CCD (this study)")
axa.set_xlim(170,0); axa.set_ylabel("CCD (m)",fontsize=PT_LABEL)
# Top at -2500 m, not -2600: the CCD peaks at -2576 m at 92-93 Ma, so the old limit
# clipped the curve at its shallowest point, which falls on OAE2.
axa.set_xlabel("Age (Ma)",fontsize=PT_LABEL); axa.set_ylim(-5100,-2500)
ax2=axa.twinx()
# The total is the series the two-predictor model uses, so it is drawn here rather than
# left implicit: heavier and neutral in colour, with the components it sums beneath it.
ax2.plot(age,D["gross_outflux"],color=COL["Total outflux"],lw=2.4,ls=(0,(5,1.6)),
         label="Total outflux (model predictor)")
ax2.plot(age,D["arc_subduction"],color=COL["Arc (subduction)"],lw=1.6,
         label="Arc outgassing (pelagic carbonates)")
ax2.plot(age,D["MOR_ridge"],color=COL["Mid-ocean ridge"],lw=1.6,
         label="Mid-ocean ridge outgassing")
ax2.plot(age,D["rift"],color=COL["Rift"],lw=1.2,label="Rift outgassing")
ax2.plot(age,D["carb_platform"],color=COL["Carbonate platform"],lw=1.2,
         label="Arc outgassing (carbonate platforms)")
# Carbon, not CO2: these fluxes are in megatonnes of CARBON per year, and the two differ
# by the molar mass ratio 44/12.
ax2.set_ylabel("Carbon outflux (Mt C yr$^{-1}$)",fontsize=PT_LABEL)
for _ax in (axa,ax2): _ax.tick_params(labelsize=PT_TICK)
# Events, ages and styling come from paper_events.py, shared with Figs 2b and 3d so the
# three panels cannot drift apart.
import importlib.util as _ilu
_espec=_ilu.spec_from_file_location("paper_events", _HERE/"paper_events.py")
_ev=_ilu.module_from_spec(_espec); _espec.loader.exec_module(_ev)
_ev.draw_events(axa, fontsize=8.5)
# Legend above the panel rather than inside it: six entries over a curve that fills the
# axes left no clear space, and a box floating on the data made the panel hard to read.
l1,la1=axa.get_legend_handles_labels(); l2,la2=ax2.get_legend_handles_labels()
# y0 = 1.07 clears the event labels, which paper_events draws just above the axes.
axa.legend(l1+l2,la1+la2,fontsize=PT_LEG,loc="lower left",
           bbox_to_anchor=(0.0,1.07,1.0,0.16),mode="expand",ncol=2,
           frameon=False,borderaxespad=0.0,handlelength=2.4,columnspacing=1.4)

# ---- (b) incremental variance ------------------------------------------------
def z(a): a=np.asarray(a,float); return (a-a.mean())/a.std()
sub=D[D["age"]<=52]; Y=z(sub["depth"])
# Sea level and CO2 outgassing only. The pelagic-carbonate sink is not a third
# independent factor: sea level is itself largely a proxy for the shift of carbonate
# burial between shelf and deep sea, so entering the sink beside it splits one
# mechanism across two predictors.
preds={"Sea level":z(sub["sl"]),"Total degassing":z(sub["gross_outflux"])}
X=np.column_stack([np.ones(len(Y))]+[preds[k] for k in preds]); b,_,_,_=np.linalg.lstsq(X,Y,rcond=None)
R2=1-np.sum((Y-X@b)**2)/np.sum((Y-Y.mean())**2)
inc={}
for drop in preds:
    keep=[k for k in preds if k!=drop]; Xk=np.column_stack([np.ones(len(Y))]+[preds[k] for k in keep])
    bk,_,_,_=np.linalg.lstsq(Xk,Y,rcond=None); inc[drop]=R2-(1-np.sum((Y-Xk@bk)**2)/np.sum((Y-Y.mean())**2))
axb1.bar(list(inc.keys()),list(inc.values()),width=0.55,
         color=["#009e73","#e69f00"],edgecolor="black")
axb1.set_ylabel("Incremental R² (unique variance)",fontsize=PT_LABEL)
axb1.set_ylim(0,0.3); axb1.set_xlim(-0.55,1.55); axb1.tick_params(labelsize=PT_TICK)
for i,(k,v) in enumerate(inc.items()):
    axb1.text(i,v+0.008,f"{v:.2f}",ha="center",va="bottom",fontsize=PT_TICK)
axb1.text(0.03,0.97,f"Two-predictor\nmodel R² = {R2:.2f}",transform=axb1.transAxes,
          ha="left",va="top",fontsize=PT_TICK)
axb1.spines["top"].set_visible(False); axb1.spines["right"].set_visible(False)

# ---- (c) correlation of each component with the CCD --------------------------
# Bars run from the lowest to the highest correlation across the model's three
# carbonate sed-rate scenarios; the marker is the mean scenario. That range is a
# sensitivity, not a confidence interval - see the note drawn in the panel.
_cc = CC.iloc[::-1].reset_index(drop=True)      # first row of the file at the top
_y = np.arange(len(_cc))
for _i, _row in _cc.iterrows():
    _c = COL.get(_row["component"], "0.4")
    axc.plot([_row["r_min_scenario"], _row["r_max_scenario"]], [_i, _i],
             color=_c, lw=3.0, solid_capstyle="butt", alpha=0.45, zorder=2)
    axc.plot([_row["r"]], [_i], marker="o", ms=7.0, zorder=3,
             mfc=_c if _row["ccd_independent"] else "white", mec=_c, mew=1.6)
axc.axvline(0.0, color="0.35", lw=0.9, zorder=1)
axc.set_yticks(_y); axc.set_yticklabels(_cc["component"], fontsize=PT_TICK)
axc.set_ylim(-0.7, len(_cc)-0.3)
axc.set_xlim(-1.0, 1.0); axc.set_xticks([-1,-0.5,0,0.5,1])
axc.set_xlabel("Correlation with CCD, 0–52 Ma", fontsize=PT_LABEL)
axc.tick_params(labelsize=PT_TICK)
axc.spines["top"].set_visible(False); axc.spines["right"].set_visible(False)
# Open markers are the components computed from the CCD itself, so their agreement with
# it is partly constructional. The p note is the substantive caveat: these series are so
# strongly autocorrelated that the effective sample size is about 3, and no correlation
# here is separable from zero.
_pmin = float(CC["p_ar1"].min())
axc.text(0.5, -0.185, f"open = computed from the CCD; all p > {np.floor(_pmin*10)/10:.1f}",
         transform=axc.transAxes, ha="center", va="top", fontsize=PT_TICK-1.0, color="0.30")

def plabel(ax,t):
    ax.text(-0.10,1.02,t,transform=ax.transAxes,fontsize=14,fontweight="bold",
            va="bottom",ha="right")
plabel(axb1,"b"); plabel(axc,"c")
# (a) sits above its legend rather than beside it, so it does not land in the legend row.
axa.text(-0.10,1.30,"a",transform=axa.transAxes,fontsize=14,fontweight="bold",
         va="bottom",ha="right")

# ---- text-collision check ----------------------------------------------------
fig.canvas.draw(); _r=fig.canvas.get_renderer()
def _live(axis,lim):
    lo,hi=min(lim),max(lim)
    return [l for loc,l in zip(axis.get_ticklocs(),axis.get_ticklabels())
            if lo<=loc<=hi and l.get_text().strip() and l.get_visible()]
_cand=list(axa.texts)+list(ax2.texts)+list(axb1.texts)+list(axc.texts)
_cand+=[axa.xaxis.label,axa.yaxis.label,ax2.yaxis.label,axb1.xaxis.label,axb1.yaxis.label,
        axc.xaxis.label,axc.yaxis.label]
_cand+=_live(axa.xaxis,axa.get_xlim())+_live(axa.yaxis,axa.get_ylim())+_live(ax2.yaxis,ax2.get_ylim())
_cand+=_live(axb1.xaxis,axb1.get_xlim())+_live(axb1.yaxis,axb1.get_ylim())
_cand+=_live(axc.xaxis,axc.get_xlim())+_live(axc.yaxis,axc.get_ylim())
_leg=axa.get_legend()
if _leg is not None: _cand+=list(_leg.get_texts())
_items=[]
for _t in _cand:
    if not (_t.get_text().strip() and _t.get_visible()): continue
    _e=_t.get_window_extent(renderer=_r)
    if _e.width>0 and _e.height>0: _items.append((_t,_e))
_bad=[]
for _i in range(len(_items)):
    for _j in range(_i+1,len(_items)):
        _a,_b=_items[_i][1],_items[_j][1]
        if _a.overlaps(_b) and min(_a.x1,_b.x1)-max(_a.x0,_b.x0)>1.0 and min(_a.y1,_b.y1)-max(_a.y0,_b.y0)>1.0:
            _bad.append(f"{_items[_i][0].get_text()!r} / {_items[_j][0].get_text()!r}")
print("  text collisions: "+(", ".join(_bad) if _bad else "none"))

for ext in ("png","pdf"): fig.savefig(f"{FIG}/Fig4_combined.{ext}",dpi=300,bbox_inches="tight")
plt.close(fig); print("wrote Fig4_combined")
