#!/usr/bin/env python3
"""Combined Figure 2 (both panels native matplotlib, identical axis width):
(a) sea-level record + RMA calibration inset (centre-bottom), Miller/Haq transition line;
(b) global CCD since 170 Ma with uncertainty + reference curves + sea level (right axis)."""
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
OUT = FIGDIR
AGEMAX=170; MHT=66.6   # Miller/Haq transition age
def rd(p,names): return pd.read_csv(p,sep=r"\s+",engine="python",comment="#",header=None,names=names)

# Sea-level curves come from step 0, i.e. normalised to GTS2020.
sl =rd(CW/"steps/step1_timescale_conversion/outputs/sealevel_shortterm_hybrid_GTS2020.txt",["age","sl"]).sort_values("age")
haq=rd(CW/"steps/step1_timescale_conversion/outputs/Haq87_longterm_GTS2020.txt",["age","sl"]).sort_values("age")
env=rd(CW/"steps/step4_sealevel_envelope/outputs/sea_level_quantile_envelope_0-205Ma.txt",["age","sl"]).sort_values("age")
mt =pd.read_csv(CW/"steps/step6_sealevel_ccd_regression/outputs/diagnostics/matched_timeseries.csv")
rr =pd.read_csv(CW/"steps/step6_sealevel_ccd_regression/outputs/diagnostics/regression_results.csv")
rma=rr[rr["method"]=="RMA"].iloc[0]
# panel b data
hyb=rd(FIGDIR/"CCD_hybrid_DM2026.txt",["age","ccd","lo","hi"]).dropna()
pred=rd(CW/"steps/step6_sealevel_ccd_regression/outputs/predicted_ccd_sl_0-205Ma.txt",["age","ccd","lo","hi"]).dropna()
bw =rd(CW/"data/reference_ccd/Boss_Wilkinson_1991_global_CCD_mean.txt",["age","ccd","mn","mx","a2","mx2","a3","mn3"]).dropna(subset=["age","ccd"])
db =rd(CW/"data/reference_ccd/Global_CCD_Delaney_Boyle.txt",["age","ccd"]).dropna(); db["ccd"]=-db["ccd"]

fig=plt.figure(figsize=(7.4,8.1))
gs=fig.add_gridspec(2,1,height_ratios=[1.0,0.98],hspace=0.20,left=0.11,right=0.865,top=0.985,bottom=0.06)

# ---------------- panel (a) ----------------
ax=fig.add_subplot(gs[0,0])
m=sl["age"]<=AGEMAX
ax.plot(sl["age"][m],sl["sl"][m],color="0.72",lw=0.7,alpha=0.9,label="Short-term record",zorder=1)
mh=haq["age"]<=AGEMAX; ax.plot(haq["age"][mh],haq["sl"][mh],color="#2b5fa6",lw=2.0,label="Haq (1987)",zorder=3)
me=env["age"]<=AGEMAX; ax.plot(env["age"][me],env["sl"][me],color="#c0392b",lw=2.4,label="Peak-following envelope",zorder=4)
ax.axhline(0,color="k",ls="--",lw=0.8,zorder=0)
ax.axvline(MHT,color="0.35",ls="--",lw=1.1,zorder=2)
# Panel geometry in mm, so the nudges below are real millimetres on the page.
_aw_mm = ax.get_position().width  * fig.get_size_inches()[0] * 25.4
_ah_mm = ax.get_position().height * fig.get_size_inches()[1] * 25.4
_dx = 1.0/_aw_mm*AGEMAX          # 1 mm along the age axis, in Myr
_dy = 2.0/_ah_mm*335.0           # 2 mm along the sea-level axis, in m
# Source of the short-term record (grey) and hence of the red envelope on each
# side of the transition. Coloured like the envelope so the link is explicit.
# Both labels are nudged 1 mm towards the transition line and 2 mm up, keeping
# them symmetric about it, so they clear the blue Haq curve without crowding each
# other across the line.
ax.text(MHT+3-_dx,262+_dy,"HaqHybrid",fontsize=9,ha="right",va="top",color="#c0392b",fontweight="bold")
ax.text(MHT-3+_dx,262+_dy,"Miller2024",fontsize=9,ha="left",va="top",color="#c0392b",fontweight="bold")
ax.set_xlim(AGEMAX,0); ax.set_ylim(-60,275)
ax.set_xlabel("Age (Ma)",fontsize=13); ax.set_ylabel("Sea level (m)",fontsize=13); ax.tick_params(labelsize=11)
ax.legend(loc="upper left",fontsize=10,frameon=False,ncol=1,handlelength=1.6,borderaxespad=0.6)
# RMA inset, centre-bottom
# Raise the inset 3 mm so its lowest tick label ("4500") clears the dashed
# zero-sea-level line, and trim 2 mm of height so the top edge still stays clear
# of the red envelope, which dips to ~123 m at 56 Ma just above the box.
iax=ax.inset_axes([0.352,0.125+3.0/_ah_mm,0.35,0.40-2.0/_ah_mm],zorder=8)
iax.set_facecolor("white"); iax.patch.set_alpha(1.0)
x=mt["Sea_level_m"].values; y=mt["CCD_m_pos_down"].values
iax.scatter(x,y,s=9,facecolors="none",edgecolors="black",linewidths=0.6,zorder=3)
xx=np.linspace(x.min(),x.max(),100); iax.plot(xx,rma["slope"]*xx+rma["intercept"],color="black",lw=1.4,zorder=4)
iax.invert_yaxis(); iax.set_xlabel("Sea level (m)",fontsize=9.5,labelpad=1.5); iax.set_ylabel("CCD (m)",fontsize=9.5,labelpad=1.5)
iax.tick_params(labelsize=8.5)
iax.text(0.04,0.96,f"CCD = {rma['intercept']:.0f} − {abs(rma['slope']):.2f}×SL\nr = {rma['r']:.2f},  R² = {rma['r2']:.2f}",
         transform=iax.transAxes,fontsize=9,va="top",ha="left",bbox=dict(boxstyle="round,pad=0.25",fc="white",ec="0.5",lw=0.5,alpha=0.8))
for s in iax.spines.values(): s.set_edgecolor("0.4")

# ---------------- panel (b) ----------------
axb=fig.add_subplot(gs[1,0])
axb.fill_between(hyb["age"],hyb["lo"],hyb["hi"],color="lightgray",zorder=1,label="CCD uncertainty")
axb.plot(pred["age"],pred["ccd"],color="cadetblue",lw=1.7,zorder=3,label="Sealevel CCD")
axb.plot(hyb["age"],hyb["ccd"],color="mediumturquoise",lw=2.6,zorder=4,label="Hybrid CCD")
axb.plot(bw["age"],bw["ccd"],color="#ee7600",lw=1.7,zorder=3,label="Boss & Wilkinson (1991)")
axb.plot(db["age"],db["ccd"],color="#8b008b",lw=1.7,zorder=3,label="Delaney & Boyle")
# paleoceanographic events discussed in the text
EVENTS=[(134.0,"WE"),(115.0,"LACI"),(93.9,"OAE2"),(56.0,"PETM"),(33.9,"EOT"),(15.2,"MMCO")]
for eage,elab in EVENTS:
    axb.axvline(eage,color="0.45",ls=":",lw=1.0,zorder=2)
    axb.text(eage,1.005,elab,transform=axb.get_xaxis_transform(),fontsize=8.5,color="0.30",
             ha="center",va="bottom",zorder=6)
axb.set_xlim(AGEMAX,0); axb.set_ylim(-5300,-2500)
axb.set_xlabel("Age (Ma)",fontsize=13); axb.set_ylabel("CCD (m)",fontsize=13); axb.tick_params(labelsize=11)
axr=axb.twinx(); axr.set_ylim(0,300)
axr.plot(env["age"],env["sl"],color="#27408b",lw=1.2,zorder=3,label="This study")
axr.plot(haq["age"],haq["sl"],color="#6ca6cd",lw=1.2,zorder=3,label="Haq (1987)")
axr.set_ylabel("Sea level (m)",fontsize=13); axr.tick_params(labelsize=11)
h1,l1=axb.get_legend_handles_labels(); h2,l2=axr.get_legend_handles_labels()
order=["Sealevel CCD","Hybrid CCD","CCD uncertainty","Boss & Wilkinson (1991)","Delaney & Boyle"]
hd={l:h for h,l in zip(h1,l1)}
leg_h=[hd[o] for o in order]+h2; leg_l=order+["Sea level: this study","Sea level: Haq (1987)"]
axb.legend(leg_h,leg_l,fontsize=7.6,frameon=False,ncol=2,loc="lower center",bbox_to_anchor=(0.5,-0.02),columnspacing=1.2,handlelength=1.5)

def plab(a,t,x=-0.085,y=1.02):
    a.text(x,y,t,transform=a.transAxes,fontsize=14,fontweight="bold",va="bottom",ha="right")
plab(ax,"a")
# panel-b label raised ~2 mm so it clears the topmost depth tick label
_bh=axb.get_position().height*fig.get_size_inches()[1]*25.4   # panel height in mm
plab(axb,"b",y=1.02+2.0/_bh)
for ext in ("png","pdf"): fig.savefig(OUT/f"Fig2_combined.{ext}",dpi=300,bbox_inches="tight")
print("wrote Fig2_combined natively")
