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
events={"OAE1a\n~120":120,"Weissert\n~133":133,"OAE2\n~93.9":93.9,"K–Pg\n66":66,
        "PETM\n56":56,"EOT\n~34":34,"MMCO\n~15":15}

fig=plt.figure(figsize=(9,10.3))
gs=fig.add_gridspec(2,2,height_ratios=[1.34,1.0],hspace=0.26,wspace=0.26)
axa=fig.add_subplot(gs[0,:]); axb1=fig.add_subplot(gs[1,0]); axb2=fig.add_subplot(gs[1,1])

# (a) CCD vs degassing
axa.plot(age,ccd,color="black",lw=2.4,label="Global CCD (this study)")
axa.set_xlim(170,0); axa.set_ylabel("CCD (m)"); axa.set_xlabel("Age (Ma)"); axa.set_ylim(-5100,-2600)
ax2=axa.twinx()
ax2.plot(age,D["arc_subduction"],color="#d1495b",lw=1.6,label="Arc outgassing (pelagic carbonates)")
ax2.plot(age,D["MOR_ridge"],color="#0072b2",lw=1.6,label="Mid-ocean ridge outgassing")
ax2.plot(age,D["rift"],color="#9467bd",lw=1.2,label="Rift outgassing")
ax2.plot(age,D["carb_platform"],color="#8c564b",lw=1.2,label="Arc outgassing (carbonate platforms)")
ax2.set_ylabel("CO$_2$ outflux (Mt C yr$^{-1}$)")
for lbl,a in events.items():
    axa.axvline(a,color="0.7",lw=0.8,ls=":"); axa.text(a,-2620,lbl,ha="center",va="top",fontsize=9.5,color="0.35")
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
preds={"Sea level":z(sub["sl"]),"Total degassing":z(sub["gross_outflux"]),"Biological sink":z(sub["sink_plate_sed"])}
X=np.column_stack([np.ones(len(Y))]+[preds[k] for k in preds]); b,_,_,_=np.linalg.lstsq(X,Y,rcond=None)
R2=1-np.sum((Y-X@b)**2)/np.sum((Y-Y.mean())**2)
inc={}
for drop in preds:
    keep=[k for k in preds if k!=drop]; Xk=np.column_stack([np.ones(len(Y))]+[preds[k] for k in keep])
    bk,_,_,_=np.linalg.lstsq(Xk,Y,rcond=None); inc[drop]=R2-(1-np.sum((Y-Xk@bk)**2)/np.sum((Y-Y.mean())**2))
axb1.bar(list(inc.keys()),list(inc.values()),color=["#009e73","#e69f00","#0072b2"],edgecolor="black")
axb1.set_ylabel("Incremental R² (unique variance)"); axb1.set_ylim(0,0.3)
axb1.tick_params(axis="x",labelrotation=15)

lags=range(-15,16); grid=age
def onlag(f,resp,mask):
    out=[]
    for L in lags:
        fs=np.interp(grid,grid+L,f,left=np.nan,right=np.nan); m=mask&np.isfinite(fs)
        out.append(np.corrcoef(fs[m],resp[m])[0,1])
    return np.array(out)
mask=(age<=52); depth=D["depth"].to_numpy()
for k,c in [("arc_subduction","#d1495b"),("sl","#009e73"),("MOR_ridge","#0072b2"),("rift","#9467bd"),("carb_platform","#8c564b")]:
    axb2.plot(list(lags),onlag(D[k].to_numpy(),depth,mask),color=c,lw=1.6,
              label={"arc_subduction":"Arc","sl":"Sea level","MOR_ridge":"MOR","gross_outflux":"Gross","rift":"Rift","carb_platform":"Carb. platform"}[k])
axb2.axvline(0,color="0.6",lw=0.8); axb2.axhline(0,color="0.6",lw=0.8)
axb2.set_xlabel("Lag (Myr; + = forcing leads CCD)"); axb2.set_ylabel("Correlation with CCD")
axb2.legend(fontsize=11,frameon=False,ncol=1,loc="center left",bbox_to_anchor=(0.0,0.52))

def plabel(ax,t):
    ax.text(-0.06,1.05,t,transform=ax.transAxes,fontsize=14,fontweight="bold",va="bottom",ha="right")
plabel(axa,"a"); plabel(axb1,"b"); plabel(axb2,"c")

for ext in ("png","pdf"): fig.savefig(f"{FIG}/Fig4_combined.{ext}",dpi=300,bbox_inches="tight")
plt.close(fig); print("wrote Fig4_combined")
