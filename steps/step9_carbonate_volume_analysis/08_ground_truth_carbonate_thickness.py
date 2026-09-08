#!/usr/bin/env python3
"""
Ground-truth the modelled present-day compacted carbonate thickness against
observed DSDP/ODP well data (Table S1), using the MEAN, MIN and MAX grids.

For each drill site with an observed total compacted carbonate thickness, the
three modelled DM2026 present-day grids (mean / min / max) are sampled at the
site, giving three modelled values and three differences (observed − modelled):

    diff_vs_mean = observed − modelled_mean   (central estimate)
    diff_vs_min  = observed − modelled_min    (largest positive; thin model)
    diff_vs_max  = observed − modelled_max    (smallest; thick model)

Outputs:
  1. Mollweide map of the modelled present-day (mean) compacted carbonate
     thickness, observed thicknesses overplotted as colour-filled circles on the
     SAME wysiwyg colour scale.                                   -> pyGMT
  2. Histogram of the central (mean) differences WITH min/max uncertainty: each
     bin shows an error bar spanning the count range obtained from the min and
     max grids.                                                  -> matplotlib
  3. Box plots of the central differences split by ocean basin
     (Pacific, Atlantic, Indian).                                -> matplotlib
  4. A CSV with observed, modelled (mean/min/max) and the three differences.

All figures are written to the shared repo Figures/ folder as 300-dpi PNG + PDF.

Default inputs (bundled in the repo under data/carbonate_thickness/):
    compacted_carbonate_thickness_present_{mean,min,max}.nc  (var 'z')
    Table_carbonate_thickness.docx
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Defaults resolve to the real project locations: the mean/min/max present-day
# grids written by the carbonate-thickness step, the benchmark table in this
# folder, and figures/CSV to output/ground_truth/.
HERE = Path(__file__).resolve().parent
SCRIPTS_ROOT = HERE.parent
FIGURES = HERE / "output" / "ground_truth"


def _thickness_root() -> Path:
    """Folder holding the carbonate_sed_thickness_*_DM2026 grid directories.

    The step was renumbered, so the current location is checked first and the
    pre-renumbering one kept as a fallback. Silently defaulting to a path that no
    longer exists aborts the carbon runner before the CO2 notebooks ever start.
    """
    for cand in (SCRIPTS_ROOT / "step8_carbonate_sediment_thickness",
                 SCRIPTS_ROOT / "step7_carbonate"):
        if (cand / "carbonate_sed_thickness_DM2026").is_dir():
            return cand
    raise SystemExit(
        "cannot find the carbonate-thickness grids: looked for "
        "carbonate_sed_thickness_DM2026/ under "
        f"{SCRIPTS_ROOT / 'step8_carbonate_sediment_thickness'} - run the "
        "carbonate-thickness step first."
    )


_GRID7 = _thickness_root()
_GRID_MEAN = _GRID7 / "carbonate_sed_thickness_DM2026" / "compacted_sediment_thickness_0.25_0.nc"
_GRID_MIN = _GRID7 / "carbonate_sed_thickness_min_DM2026" / "compacted_sediment_thickness_0.25_0.nc"
_GRID_MAX = _GRID7 / "carbonate_sed_thickness_max_DM2026" / "compacted_sediment_thickness_0.25_0.nc"
_TABLE = HERE / "Table_carbonate_thickness.docx"


def _paper_figures_dir():
    """The paper's figure folder, so Figs S2-S4 land beside Figs 1-4 and S1.

    This script is run standalone by the carbon runner, from its own directory, so
    the workflow package may not be importable; fall back to locating a sibling
    Paper/Figures by walking up from here.
    """
    try:
        from ccdworkflow import config           # noqa: PLC0415
        if config.PAPER_FIGURES is not None:
            return config.PAPER_FIGURES
    except Exception:
        pass
    for base in (HERE, *HERE.parents):
        cand = base.parent / "Paper" / "Figures"
        if cand.is_dir():
            return cand
    return None


PAPER_FIGURES = _paper_figures_dir()


def _targets(name, paper_name):
    """Where a figure is written: this step's output folder, plus the paper's figure
    folder under its supplement number when the figure is one of the supplement's."""
    FIGURES.mkdir(parents=True, exist_ok=True)
    out = [(FIGURES, name)]
    if paper_name:
        if PAPER_FIGURES is None:
            print(f"  [figures] paper figure folder not found; {paper_name} not propagated")
        else:
            PAPER_FIGURES.mkdir(parents=True, exist_ok=True)
            out.append((PAPER_FIGURES, paper_name))
    return out


def _save_mpl(fig, name, dpi=300, paper_name=None):
    for folder, stem in _targets(name, paper_name):
        fig.savefig(folder / f"{stem}.png", dpi=dpi, bbox_inches="tight")
        fig.savefig(folder / f"{stem}.pdf", dpi=dpi, bbox_inches="tight")


def _save_pygmt(fig, name, dpi=300, paper_name=None):
    for folder, stem in _targets(name, paper_name):
        fig.savefig(str(folder / f"{stem}.pdf"))
        fig.savefig(str(folder / f"{stem}.png"), dpi=dpi)


# --------------------------------------------------------------------------
# Observed table (Table S1)
# --------------------------------------------------------------------------
def read_table(docx_path: Path) -> pd.DataFrame:
    import docx
    tbl = docx.Document(str(docx_path)).tables[0]
    rows = []
    for r in tbl.rows[2:]:                       # first two rows are the header
        c = [x.text.strip().replace("\n", " ") for x in r.cells]
        try:
            rows.append({"Site": c[0], "lon": float(c[1]), "lat": float(c[2]),
                         "region": c[3], "observed_m": float(c[8])})
        except (ValueError, IndexError):
            continue
    return pd.DataFrame(rows)


def basin3(lon: float) -> str:
    """Coarse 3-way basin from longitude: Atlantic / Indian / Pacific."""
    if -70 <= lon < 20:
        return "Atlantic"
    if 20 <= lon < 145:
        return "Indian"
    return "Pacific"


# --------------------------------------------------------------------------
# Sample a grid at each site (bilinear; nearest-finite fallback)
# --------------------------------------------------------------------------
def sample_grid(da: xr.DataArray, lons, lats, search_deg: float = 1.5) -> np.ndarray:
    latv, lonv = da["lat"].values, da["lon"].values
    out = np.full(len(lons), np.nan)
    interp = da.interp(lat=("pts", np.asarray(lats)), lon=("pts", np.asarray(lons)),
                       method="linear").values
    for i, (lo, la, v) in enumerate(zip(lons, lats, interp)):
        if np.isfinite(v):
            out[i] = v
            continue
        jl = np.where(np.abs(latv - la) <= search_deg)[0]
        jo = np.where(np.abs(((lonv - lo + 180) % 360) - 180) <= search_deg)[0]
        if jl.size and jo.size:
            sub = da.values[np.ix_(jl, jo)]
            if np.isfinite(sub).any():
                la_g, lo_g = np.meshgrid(latv[jl], lonv[jo], indexing="ij")
                d = (la_g - la) ** 2 + (((lo_g - lo + 180) % 360 - 180)) ** 2
                d[~np.isfinite(sub)] = np.inf
                out[i] = sub.flat[np.argmin(d)]
    return out


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def make_map(da_mean, df, vmax, name, paper_name=None):
    import pygmt
    tmp = Path("/tmp/_gt_grid.nc")
    da_mean.to_netcdf(tmp)
    pygmt.config(FONT_ANNOT_PRIMARY="10p,Helvetica", FONT_LABEL="12p,Helvetica")
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="wysiwyg", series=[0, vmax], background=True)
    fig.grdimage(grid=str(tmp), region="d", projection="W15c", cmap=True,
                 frame=["xa60f30", "ya30f15", "WSne"])
    fig.coast(shorelines="0.3p,gray30", land="gray85", region="d", projection="W15c")
    # Double outline for visibility on any background colour: plot three stacked
    # circles of decreasing size — white halo (outer), black ring, then the
    # value-filled circle on top.
    fig.plot(x=df["lon"], y=df["lat"], fill="white", style="c0.32c", pen="0.25p,white")
    fig.plot(x=df["lon"], y=df["lat"], fill="black", style="c0.28c", pen="0.25p,black")
    fig.plot(x=df["lon"], y=df["lat"], fill=df["observed_m"], cmap=True, style="c0.24c")
    ann = max(round((vmax / 4) / 10) * 10, 10)
    fig.colorbar(frame=f"xa{ann}f{ann/2:g}+lCompacted carbonate thickness (m)",
                 position="JBC+w10c/0.4c+o0/1c")
    _save_pygmt(fig, name, dpi=300, paper_name=paper_name)
    tmp.unlink(missing_ok=True)


def make_histogram(df, name, paper_name=None):
    """Histogram of central (mean) differences with min/max grid uncertainty.

    Bars = site counts per bin using the MEAN grid. Error bars span the count
    range obtained from the MIN and MAX grids (the two ends of the modelled
    carbonate-thickness envelope), so each bin shows how many sites could move
    in/out of it across the model uncertainty.
    """
    bins = np.arange(-90, 91, 10)
    centres = (bins[:-1] + bins[1:]) / 2.0
    width = (bins[1] - bins[0]) * 0.9

    c_mean, _ = np.histogram(df["diff_vs_mean"].dropna(), bins=bins)
    c_min, _ = np.histogram(df["diff_vs_min"].dropna(), bins=bins)
    c_max, _ = np.histogram(df["diff_vs_max"].dropna(), bins=bins)

    d = df["diff_vs_mean"].dropna()
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    # Central estimate: filled bars from the mean grid.
    ax.bar(centres, c_mean, width=width, color="#4C72B0", edgecolor="black",
           alpha=0.85, label="Central (mean grid)")
    # Uncertainty: the min- and max-grid distributions as step outlines, so bins
    # that gain/lose sites under the thinner/thicker grid are shown as the
    # envelope of the distribution rather than as detached error bars.
    ax.stairs(c_min, bins, color="#D55E00", lw=1.6, label="min grid (thinner)")
    ax.stairs(c_max, bins, color="#009E73", lw=1.6, label="max grid (thicker)")
    ax.axvline(0, color="black", lw=1)
    ax.axvline(d.mean(), color="crimson", lw=1.6, ls="--", label=f"mean = {d.mean():+.1f} m")
    ax.axvline(d.median(), color="darkgreen", lw=1.6, ls=":", label=f"median = {d.median():+.1f} m")
    ax.set_xlabel("Observed − modelled compacted carbonate thickness (m)")
    ax.set_ylabel("Number of sites")
    # No title: this is Fig. S3 of the supplement, and n and RMS are given in its caption.
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    _save_mpl(fig, name, dpi=300, paper_name=paper_name)
    plt.close(fig)


def make_boxplot(df, name, paper_name=None):
    order = ["Global", "Pacific", "Atlantic", "Indian"]
    colors = {"Global": "#7F7F7F", "Pacific": "#D55E00", "Atlantic": "#0072B2", "Indian": "#009E73"}
    data = [df["diff_vs_mean"].dropna().values if b == "Global"
            else df.loc[df["basin"] == b, "diff_vs_mean"].dropna().values
            for b in order]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bp = ax.boxplot(data, tick_labels=[f"{b}\n(n={len(d)})" for b, d in zip(order, data)],
                    patch_artist=True, showmeans=True, widths=0.6,
                    medianprops=dict(color="black", lw=1.5),
                    meanprops=dict(marker="D", markerfacecolor="white", markeredgecolor="black"))
    for patch, b in zip(bp["boxes"], order):
        patch.set_facecolor(colors[b]); patch.set_alpha(0.7)
    for i, d in enumerate(data, 1):
        ax.scatter(np.random.default_rng(i).normal(i, 0.05, len(d)), d,
                   s=18, color="black", alpha=0.5, zorder=3)
    ax.axhline(0, color="black", lw=1)
    ax.set_ylabel("Observed − modelled thickness (m)")
    # No title: this is Fig. S4 of the supplement; the caption describes it.
    fig.tight_layout()
    _save_mpl(fig, name, dpi=300, paper_name=paper_name)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=str(_TABLE))
    ap.add_argument("--grid-mean", default=str(_GRID_MEAN))
    ap.add_argument("--grid-min", default=str(_GRID_MIN))
    ap.add_argument("--grid-max", default=str(_GRID_MAX))
    ap.add_argument("--outdir", default=str(HERE / "output" / "ground_truth"))
    ap.add_argument("--vmax", type=float, default=150.0,
                    help="shared colour-scale max (m); matched to the well-data range")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    df = read_table(Path(args.table))

    grids = {k: xr.open_dataset(p)["z"]
             for k, p in [("mean", args.grid_mean), ("min", args.grid_min), ("max", args.grid_max)]}
    for k, da in grids.items():
        df[f"modelled_{k}_m"] = sample_grid(da, df["lon"].values, df["lat"].values)

    df["diff_vs_mean"] = df["observed_m"] - df["modelled_mean_m"]
    df["diff_vs_min"] = df["observed_m"] - df["modelled_min_m"]
    df["diff_vs_max"] = df["observed_m"] - df["modelled_max_m"]
    df["basin"] = df["lon"].map(basin3)

    csv = outdir / "carbonate_thickness_observed_vs_modelled.csv"
    df.to_csv(csv, index=False, float_format="%.2f")

    # The paper_name arguments send Figs S2-S4 straight into the paper's figure
    # folder alongside Figs 1-4 and S1, so the supplement can be assembled from one place.
    make_map(grids["mean"], df, args.vmax, "map_obs_vs_modelled_carbonate_thickness",
             paper_name="FigS2_modelled_vs_observed_thickness")
    make_histogram(df, "difference_histogram",
                   paper_name="FigS3_thickness_residual_histogram")
    make_boxplot(df, "difference_boxplot_by_ocean",
                 paper_name="FigS4_thickness_residuals_by_basin")

    print(f"n sites: {len(df)}")
    for col, lbl in [("diff_vs_mean", "mean grid"), ("diff_vs_min", "min grid"), ("diff_vs_max", "max grid")]:
        s = df[col].dropna()
        print(f"  obs-modelled ({lbl:9s}): mean {s.mean():+6.1f}  median {s.median():+6.1f}  RMS {np.sqrt((s**2).mean()):5.1f} m")
    print("  by basin (central/mean grid):")
    for b in ["Pacific", "Atlantic", "Indian"]:
        s = df.loc[df.basin == b, "diff_vs_mean"].dropna()
        print(f"    {b:9s} n={len(s):2d}  mean {s.mean():+6.1f}  median {s.median():+6.1f} m")
    print(f"wrote: {csv}  and figures -> {FIGURES}")


if __name__ == "__main__":
    main()
