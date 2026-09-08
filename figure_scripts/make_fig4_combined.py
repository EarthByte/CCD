#!/usr/bin/env python3
"""Fig. 4 (combined, vertically stacked):
   (a) Global CCD vs recomputed solid-Earth CO2 outflux components (arc, MOR, gross);
   (b) attribution — incremental variance (left) and lag correlation (right), 0-52 Ma.
Panels carry boxed lower-case labels. Reads attribution_matched_series.csv next to this script."""
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
age=D["age"].to_numpy(); ccd=-D["depth"].to_numpy()

# Two panels. The lag panel is gone: carbonate compensation restores ocean saturation
# in about 5,000-10,000 years, whether the perturbation is a sea-level-driven shift in
# carbonate burial or an added solid-Earth CO2 input, which is 100 to 200 times shorter
# than the 1 Myr sampling. A lead or lag of several Myr is not resolvable here, so only
# the zero-lag relationship is shown.
fig=plt.figure(figsize=(7.48,8.6))
gs=fig.add_gridspec(2,1,height_ratios=[1.34,0.95],hspace=0.30,top=0.885)
axa=fig.add_subplot(gs[0,0]); axb1=fig.add_subplot(gs[1,0])

# (a) CCD vs degassing
axa.plot(age,ccd,color="black",lw=2.4,label="Global CCD (this study)")
axa.set_xlim(170,0); axa.set_ylabel("CCD (m)"); axa.set_xlabel("Age (Ma)"); axa.set_ylim(-5100,-2600)
ax2=axa.twinx()
ax2.plot(age,D["arc_subduction"],color="#d1495b",lw=1.6,label="Arc outgassing (pelagic carbonates)")
ax2.plot(age,D["MOR_ridge"],color="#0072b2",lw=1.6,label="Mid-ocean ridge outgassing")
ax2.plot(age,D["rift"],color="#9467bd",lw=1.2,label="Rift outgassing")
ax2.plot(age,D["carb_platform"],color="#8c564b",lw=1.2,label="Arc outgassing (carbonate platforms)")
ax2.set_ylabel("CO$_2$ outflux (Mt C yr$^{-1}$)")
# Events, ages and styling come from paper_events.py, shared with Figs 2b and 3d so the
# three panels cannot drift apart. They had: this panel previously carried a different
# set, with the Valanginian event at 133 Ma against 134 Ma in Fig. 2b, an early Aptian
# event at 120 Ma against 115 Ma, and a K-Pg line the other panels lacked.
import importlib.util as _ilu
_espec=_ilu.spec_from_file_location("paper_events", _HERE/"paper_events.py")
_ev=_ilu.module_from_spec(_espec); _espec.loader.exec_module(_ev)
_ev.draw_events(axa, fontsize=8.5)
l1,la1=axa.get_legend_handles_labels(); l2,la2=ax2.get_legend_handles_labels()
# Legend on the LEFT, with the top of the box at -3450 m on the CCD axis. Set by
# depth rather than by axes fraction so it stays put if the y-limits change.
_LEG_TOP_M = -3450.0
_y0,_y1 = axa.get_ylim()
axa.legend(l1+l2,la1+la2,fontsize=10,loc="upper left",
           bbox_to_anchor=(0.012,(_LEG_TOP_M-_y0)/(_y1-_y0)),
           ncol=1,frameon=True,facecolor="white",framealpha=0.72,
           edgecolor="0.5").set_zorder(20)

# (b) attribution
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
axb1.set_ylabel("Incremental R² (unique variance)")
axb1.set_ylim(0,0.3); axb1.set_xlim(-0.55,1.55)
for i,(k,v) in enumerate(inc.items()):
    axb1.text(i,v+0.008,f"{v:.2f}",ha="center",va="bottom",fontsize=10)
axb1.text(0.98,0.94,f"Two-predictor model R² = {R2:.2f}",transform=axb1.transAxes,
          ha="right",va="top",fontsize=10)
axb1.spines["top"].set_visible(False); axb1.spines["right"].set_visible(False)

def plabel(ax,t):
    ax.text(-0.06,1.05,t,transform=ax.transAxes,fontsize=14,fontweight="bold",va="bottom",ha="right")
plabel(axa,"a"); plabel(axb1,"b")

fig.canvas.draw(); _r=fig.canvas.get_renderer()
def _live(axis,lim):
    lo,hi=min(lim),max(lim)
    return [l for loc,l in zip(axis.get_ticklocs(),axis.get_ticklabels())
            if lo<=loc<=hi and l.get_text().strip() and l.get_visible()]
_cand=list(axa.texts)+list(ax2.texts)+list(axb1.texts)
_cand+=[axa.xaxis.label,axa.yaxis.label,ax2.yaxis.label,axb1.xaxis.label,axb1.yaxis.label]
_cand+=_live(axa.xaxis,axa.get_xlim())+_live(axa.yaxis,axa.get_ylim())+_live(ax2.yaxis,ax2.get_ylim())
_cand+=_live(axb1.xaxis,axb1.get_xlim())+_live(axb1.yaxis,axb1.get_ylim())
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
