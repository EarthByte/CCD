#!/usr/bin/env python3
"""
Attribution (idea 1) + carbonate-budget closure (idea 2) analysis.

Inputs (staged):
  CCD hybrid curve      : <figures>/CCD_hybrid_DM2026.txt  (age, CCD, min, max; CCD negative-up)
  Sea level (long-term) : step3 sea_level_quantile_envelope_0-205Ma.txt
  Degassing components  : 05_atmospheric_influx_all_sources.csv (Mt C/yr, per component, min/mean/max)
  Plate influx by reservoir : 02_plate_influx.csv (Mt C/yr; used for the biological-sink proxy)

Conventions:
  depth = -CCD  (positive down; larger = deeper = higher deep-ocean saturation)
  Primary independent window = 0-52 Ma (CCD is data-driven there; degassing is
  independent of how the CCD was built). 0-170 Ma reported for context only.
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")

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
SL  = CW / "steps/step4_sealevel_envelope/outputs/sea_level_quantile_envelope_0-205Ma.txt"
# degassing components from the consolidated headless CO2 notebooks (Notebook 05):
ATM = CW / "steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Outputs/Notebook05/csv/05_atmospheric_influx_all_sources.csv"
# plate influx split by reservoir (Notebook 02), needed for the biological-sink proxy
PLI = CW / "steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Outputs/Notebook02/csv/02_plate_influx.csv"

# ---------- load ----------
ccd = pd.read_csv(FIGDIR / "CCD_hybrid_DM2026.txt", sep=r"\s+", header=None,
                  names=["age","ccd","cmin","cmax"]).dropna()
ccd = ccd[ccd["age"]<=170].copy()
ccd["depth"] = -ccd["ccd"]                      # positive down

sl = pd.read_csv(SL, sep=r"\s+|\t+", engine="python", comment="#", header=None, names=["age","sl"]).dropna()

atm = pd.read_csv(ATM)
atm = atm.rename(columns={"Age (Ma)":"age"})

grid = np.arange(0,171,1.0)
def onto(a, v):
    a=np.asarray(a,float); v=np.asarray(v,float); o=np.argsort(a); a,v=a[o],v[o]
    keep=np.concatenate([[True],np.diff(a)>0]); return np.interp(grid,a[keep],v[keep])

D = pd.DataFrame({"age":grid})
D["depth"]  = onto(ccd["age"], ccd["depth"])
D["sl"]     = onto(sl["age"], sl["sl"])
comp = {
 "MOR_ridge":       "ridge_outflux_mean",
 "arc_subduction":  "subduction_outflux_mean",
 "rift":            "rift_outflux_biased_mean",
 "carb_platform":   "carbonate_platform_outflux_mean",
 "intraplate":      "intraplate_volcanism_outflux_mean",
 "gross_outflux":   "gross_atmospheric_outflux_biased_rift_mean",
 "net_influx":      "net_atmospheric_influx_biased_and_sed_mean",
 # the two end-members of the net solid-Earth outflux (Fig. S5 and Table S1):
 # counting the growing pelagic carbonate reservoir as a genuinely new sink
 # (net_out_incl_sed) or as a shallow-to-deep redistribution of burial
 # (net_out_excl_sed, read separately below).
 "net_out_incl_sed": "net_atmospheric_influx_unbiased_and_sed_mean",
}
for k,c in comp.items():
    D[k]=onto(atm["age"], atm[c]) if c in atm.columns else np.nan

NET = CW / "steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Outputs/Notebook05/csv/05_net_carbon_outflux.csv"
_net = pd.read_csv(NET, index_col=0, header=[0, 1])
_net.index = pd.to_numeric(_net.index, errors="coerce")
_net = _net[np.isfinite(_net.index)]
D["net_out_excl_sed"] = onto(_net.index.to_numpy(float), _net[("net_carbon_outflux", "mean")].to_numpy())

# ---------- biological-sink proxy ----------
# Carbon carried down on the subducting plate and retained by the upper plate.
# The model's `gross_upper_plate_influx_with_sed` sums sediments + altered crust +
# SERPENTINITE, but serpentinite carbon sits in serpentinised mantle lithosphere:
# it is mantle-derived and returns to the mantle, a closed mantle-plate-mantle
# loop, so it is not a biological sink and is excluded here. Altered oceanic
# crust is kept: its carbon is hydrothermal, but ultimately of biological origin.
# Mantle lithosphere proper was never in the sum. Serpentinite is ~5% of the old
# total over 0-52 Ma, so this is a small correction of principle.
_pli = pd.read_csv(PLI, index_col=0, header=[0, 1])
_pli.index = pd.to_numeric(_pli.index, errors="coerce")
_pli = _pli[np.isfinite(_pli.index)]
_bio = _pli[("sediments", "mean")].to_numpy() + _pli[("crust", "mean")].to_numpy()
D["sink_plate_sed"] = onto(_pli.index.to_numpy(float), _bio)

# ---------- AR1-adjusted correlation ----------
def ar1_corr(x,y):
    m=np.isfinite(x)&np.isfinite(y); x,y=x[m],y[m]; n=len(x)
    if n<6: return np.nan,np.nan,np.nan
    r,_=stats.pearsonr(x,y)
    r1x=stats.pearsonr(x[:-1],x[1:]).statistic; r1y=stats.pearsonr(y[:-1],y[1:]).statistic
    neff=float(np.clip(n*(1-r1x*r1y)/(1+r1x*r1y),3,n))
    t=r*np.sqrt((neff-2)/max(1e-12,1-r*r)); p=2*stats.t.sf(abs(t),df=neff-2)
    return float(r),neff,float(p)

# No lag scan. Carbonate compensation restores the ocean's saturation state in roughly
# 5,000-10,000 years, whether the perturbation is a sea-level-driven shift in carbonate
# burial or an added solid-Earth CO2 input. That is 100 to 200 times shorter than the
# 1 Myr sampling interval, so on this grid the CCD responds within a single step and a
# lead or lag of several Myr is not a physical quantity these series can resolve. Only
# the zero-lag correlation is reported.

for label, mask in [("OBSERVED 0-52 Ma", (D["age"]<=52).to_numpy()),
                    ("FULL 0-170 Ma (context; >52 model)", (D["age"]<=170).to_numpy())]:
    depthF=D["depth"].to_numpy(); ddepthF=np.gradient(depthF)
    depth=depthF[mask]; ddepth=ddepthF[mask]
    print(f"\n================ {label}  (n={int(mask.sum())}) ================")
    print(f"{'series':16s} {'r(CCD)':>11s} {'pAR1':>7s} | {'r(dCCD)':>8s}")
    for k in ["sl","MOR_ridge","arc_subduction","rift","carb_platform","intraplate",
              "gross_outflux","net_influx","net_out_excl_sed","net_out_incl_sed","sink_plate_sed"]:
        sF=D[k].to_numpy()
        r0,_,p0=ar1_corr(sF[mask],depth)
        rd,_,_=ar1_corr(np.gradient(sF)[mask],ddepth)
        print(f"{k:16s} {r0:11.2f} {p0:7.3f} | {rd:8.2f}")

# ---------- multiple regression attribution (0-52, standardized) ----------
sub=D[D["age"]<=52].copy()
def z(a): a=np.asarray(a,float); return (a-a.mean())/a.std()
Y=z(sub["depth"])
# Two predictors only. The pelagic-carbonate sink is not a third independent factor:
# sea level is itself largely a proxy for the shift of carbonate burial between shelf
# and deep sea, so entering the sink alongside it splits one mechanism in two. It is
# still reported as a correlation above.
preds={"sea_level":z(sub["sl"]), "total_degassing":z(sub["gross_outflux"])}
X=np.column_stack([preds[k] for k in preds]); X=np.column_stack([np.ones(len(Y)),X])
beta,_,_,_=np.linalg.lstsq(X,Y,rcond=None)
yhat=X@beta; R2=1-np.sum((Y-yhat)**2)/np.sum((Y-Y.mean())**2)
print("\n================ MULTIPLE REGRESSION (0-52 Ma, standardized) ================")
print(f"CCD ~ sea_level + total_degassing   R2_full={R2:.2f}")
for name,b in zip(["intercept"]+list(preds),beta):
    print(f"  {name:16s} std beta = {b:+.2f}")
# incremental R2 (drop-one)
for drop in preds:
    keep=[k for k in preds if k!=drop]
    Xk=np.column_stack([np.ones(len(Y))]+[preds[k] for k in keep])
    bk,_,_,_=np.linalg.lstsq(Xk,Y,rcond=None); yk=Xk@bk
    R2k=1-np.sum((Y-yk)**2)/np.sum((Y-Y.mean())**2)
    print(f"  incremental R2 of {drop:16s} = {R2-R2k:+.2f}")

# ---------- per-component correlations, with the sed-rate scenario spread ----------
# Table S1 and Figure 4c both read this file, so the numbers cannot drift apart.
# The spread is the model's own min/mean/max carbonate sed-rate scenarios. It is a
# sensitivity range, NOT a confidence interval: the statistical interval is useless
# here, because the CCD and the fluxes are so strongly autocorrelated that the
# effective sample size falls to its floor and the 95% interval on every one of these
# correlations covers -1 to +1. That is the substantive result, and is why none of
# them is significant once the autocorrelation is accounted for.
_COMPONENTS = [
    ("Mid-ocean ridge",    "ridge_outflux",                       True),
    ("Rift",               "rift_outflux_biased",                 True),
    ("Carbonate platform", "carbonate_platform_outflux",          True),
    ("Arc (subduction)",   "subduction_outflux",                  False),
    ("Total outflux",      "gross_atmospheric_outflux_biased_rift", False),
]
_rows = []
_m52 = (D["age"] <= 52).to_numpy()
_depth52 = D["depth"].to_numpy()[_m52]
for _label, _base, _indep in _COMPONENTS:
    _r = {}
    for _sc in ("mean", "min", "max"):
        _col = f"{_base}_{_sc}"
        if _col not in atm.columns:
            _r[_sc] = (np.nan, np.nan, np.nan); continue
        _y = onto(atm["age"], atm[_col])[_m52]
        _r[_sc] = ar1_corr(_depth52, _y)
    # Range across all three scenarios, the mean included: a min or max sed-rate grid
    # does not have to bracket the mean one, and on the figure the marker has to sit
    # inside its own bar.
    _all = [_r[_sc][0] for _sc in ("mean", "min", "max")]
    _rows.append({"component": _label, "ccd_independent": _indep,
                  "r": _r["mean"][0], "r_min_scenario": min(_all), "r_max_scenario": max(_all),
                  "n_effective": _r["mean"][1], "p_ar1": _r["mean"][2]})
_cc = pd.DataFrame(_rows)
_cc.to_csv(_HERE / "component_correlations.csv", index=False)
print("\n================ COMPONENT CORRELATIONS (0-52 Ma) ================")
for _row in _rows:
    print(f"  {_row['component']:20s} r = {_row['r']:+.2f}  "
          f"scenarios [{_row['r_min_scenario']:+.2f}, {_row['r_max_scenario']:+.2f}]  "
          f"n_eff = {_row['n_effective']:.1f}  p = {_row['p_ar1']:.2f}"
          f"{'' if _row['ccd_independent'] else '   (CCD-derived)'}")

D.to_csv(_HERE / "attribution_matched_series.csv", index=False)
print("\nsaved matched series ->", _HERE / "attribution_matched_series.csv")
