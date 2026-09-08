#!/usr/bin/env python3
"""
Figure S5 - distribution of carbonate accumulation rate with water depth, at
80, 60, 40, 20 and 0 Ma.

What is plotted
---------------
For every abyssal grid cell at reconstruction time T:

    accumulation rate = compacted carbonate thickness(T) / seafloor age(T)

in m/Myr, i.e. the mean rate at which carbonate has accumulated on that crust over
its lifetime up to T. Cells are grouped into 500 m water-depth bins using the
paleobathymetry at the same time, and the distribution within each bin is drawn as
a violin. The CCD at T is marked, so the figure shows accumulation collapsing at
the compensation depth and that collapse migrating downward as the CCD deepens.

Why not a difference of two thickness grids
-------------------------------------------
The grids are reconstructed into palaeo-coordinates, so a given cell holds
different crust at different times: a cell-by-cell difference between the grids at
T and T+20 subtracts unrelated pieces of seafloor. Dividing a single time slice by
seafloor age needs no differencing and stays within one reconstruction.

The rate is therefore a lifetime average, not an instantaneous one, and it is
plotted against the depth the crust occupies at T rather than the depth at which
each increment of carbonate was laid down. Both are stated in the caption.

Inputs (all already in the workflow):
    steps/step8_carbonate_sediment_thickness/carbonate_sed_thickness_DM2026/
        compacted_sediment_thickness_0.25_{T}.nc
    steps/step8_carbonate_sediment_thickness/input_grids/
        Alfonso2024_pybacktrack_merged_paleobathymetry/paleobathymetry_{T}Ma.nc
        Alfonso2024_SeafloorAgeGrids-2026/SEAFLOOR_AGE_grid_{T:.2f}Ma.nc
    Figures/CCD_hybrid_DM2026.txt                     (the CCD at each age)

Output: Figures/FigS5_accumulation_rate_vs_depth.{png,pdf}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- path resolution (same convention as the other Paper/ scripts) ----------
_HERE = Path(__file__).resolve().parent


def _workflow_root() -> Path:
    for _p in (_HERE.parent / "CCD_workflow_clean", _HERE.parent, _HERE):
        if (_p / "steps").is_dir() and (_p / "ccdworkflow").is_dir():
            return _p
    raise SystemExit("cannot locate the CCD workflow root (expected steps/ + ccdworkflow/)")


CW = _workflow_root()
FIGDIR = _HERE / "Figures" if (_HERE / "Figures").is_dir() else CW / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)
STEP8 = CW / "steps/step8_carbonate_sediment_thickness"
CARB = STEP8 / "carbonate_sed_thickness_DM2026"
GRIDS = STEP8 / "input_grids"
BATHY = GRIDS / "Alfonso2024_pybacktrack_merged_paleobathymetry"
AGES = GRIDS / "Alfonso2024_SeafloorAgeGrids-2026"
HYBRID = FIGDIR / "CCD_hybrid_DM2026.txt"

TIMES = (80, 60, 40, 20, 0)
DEPTH_EDGES = np.arange(1000.0, 6001.0, 500.0)     # m below sea level
MIN_AGE_MYR = 1.0          # skip ridge-crest cells, where thickness/age is unstable
MIN_CELLS = 400            # a violin needs enough cells to mean anything
MAX_PER_BIN = 20000        # subsample for the kernel density estimate
RATE_MAX = 8.0             # m/Myr, axis limit
SEED = 12345

# Geology artwork: nothing below ~7 pt at final size.
PT_TICK, PT_LABEL, PT_PANEL, PT_ANNOT = 7.5, 8.5, 10.0, 7.5


def _to_pm180(da: xr.DataArray) -> xr.DataArray:
    """Grids come in both 0..360 and -180..180; normalise, dropping the wrap column."""
    if float(da["lon"].max()) > 180.0:
        da = da.assign_coords(lon=((da.lon + 180.0) % 360.0) - 180.0).sortby("lon")
        _, keep = np.unique(da["lon"].values, return_index=True)
        da = da.isel(lon=np.sort(keep))
    return da


def _open(path: Path, var: str = "z") -> xr.DataArray:
    if not path.exists():
        raise SystemExit(f"missing input: {path}")
    return xr.open_dataset(path)[var]


def sample(T: int):
    """Depth (m, positive down) and accumulation rate (m/Myr) for every abyssal cell."""
    carb = _open(CARB / f"compacted_sediment_thickness_0.25_{T}.nc")
    bath = _to_pm180(_open(BATHY / f"paleobathymetry_{T}Ma.nc"))
    age = _to_pm180(_open(AGES / f"SEAFLOOR_AGE_grid_{T:.2f}Ma.nc"))
    B = bath.interp(lon=carb.lon, lat=carb.lat, method="linear").values
    A = age.interp(lon=carb.lon, lat=carb.lat, method="nearest").values
    C = carb.values
    m = np.isfinite(C) & np.isfinite(B) & np.isfinite(A) & (B < 0.0) & (A > MIN_AGE_MYR)
    return -B[m], C[m] / A[m]


def ccd_curve():
    h = pd.read_csv(HYBRID, sep=r"\s+", header=None, names=["age", "ccd", "lo", "hi"])
    return lambda t: -float(np.interp(t, h["age"], h["ccd"]))     # positive down


def draw_violin(ax, sel, centre, colour, width=400.0, npts=256):
    """A horizontal violin whose density is evaluated only over the data range.

    matplotlib's violinplot pads the kernel support beyond the observed values, which
    on these distributions draws a long hairline tail across the whole panel and, at
    the deep bins where almost every cell is zero, extends the density below zero.
    Rates cannot be negative, so the support is clipped to [min, max] of the data.
    """
    lo, hi = float(sel.min()), float(sel.max())
    if hi <= lo:
        return
    kde = gaussian_kde(sel)
    x = np.linspace(lo, hi, npts)
    dens = kde(x)
    dens = dens / dens.max() * (width / 2.0)
    ax.fill_betweenx(centre + np.concatenate([dens, -dens[::-1]]),
                     np.concatenate([x, x[::-1]]),
                     facecolor=colour, edgecolor="0.25", linewidth=0.4, zorder=3)
    med = float(np.median(sel))
    d_med = float(np.interp(med, x, dens))
    ax.plot([med, med], [centre - d_med, centre + d_med], color="black", lw=1.0, zorder=4)


def check_text_collisions(fig, axes) -> list[str]:
    """Report overlapping text, rather than relying on anyone spotting it by eye."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = []
    for ax in axes:
        for t in ax.texts + ax.get_xticklabels() + ax.get_yticklabels() + [ax.xaxis.label, ax.yaxis.label]:
            if t.get_text().strip():
                items.append((t, t.get_window_extent(renderer=r)))
    bad = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][1], items[j][1]
            if a.overlaps(b):
                ov = min(a.x1, b.x1) - max(a.x0, b.x0), min(a.y1, b.y1) - max(a.y0, b.y0)
                if ov[0] > 1.0 and ov[1] > 1.0:
                    bad.append(f"{items[i][0].get_text()!r} overlaps {items[j][0].get_text()!r}")
    return bad


def main() -> None:
    rng = np.random.default_rng(SEED)
    ccd_at = ccd_curve()
    centres = 0.5 * (DEPTH_EDGES[:-1] + DEPTH_EDGES[1:])
    # sequential single hue: deeper water = darker blue
    shades = plt.get_cmap("Blues")(np.linspace(0.32, 0.92, len(centres)))

    fig, axes = plt.subplots(1, len(TIMES), figsize=(7.48, 4.3), sharey=True)
    for k, (ax, T) in enumerate(zip(axes, TIMES)):
        depth, rate = sample(T)
        n_tot = 0
        for i, (lo, hi) in enumerate(zip(DEPTH_EDGES[:-1], DEPTH_EDGES[1:])):
            sel = rate[(depth >= lo) & (depth < hi)]
            if sel.size < MIN_CELLS:
                continue
            n_tot += sel.size
            if sel.size > MAX_PER_BIN:
                sel = rng.choice(sel, MAX_PER_BIN, replace=False)
            if np.allclose(sel, sel[0]):            # a constant bin has no density
                continue
            draw_violin(ax, sel, centres[i], shades[i])

        c = ccd_at(T)
        ax.axhline(c, color="#b3202c", lw=1.2, ls="--", zorder=5)
        ax.text(RATE_MAX * 0.97, c - 90, "CCD", color="#b3202c", fontsize=PT_ANNOT,
                ha="right", va="bottom")
        # Panel letter and age sit above the frame, so neither can ever land on a violin.
        ax.text(-0.02, 1.02, "abcde"[k], transform=ax.transAxes, fontsize=PT_PANEL,
                fontweight="bold", ha="left", va="bottom")
        ax.text(1.0, 1.03, f"{T} Ma", transform=ax.transAxes, fontsize=PT_ANNOT,
                ha="right", va="bottom")
        ax.set_xlim(0, RATE_MAX)
        ax.set_ylim(DEPTH_EDGES[-1], DEPTH_EDGES[0])
        ax.set_xticks([0, 2, 4, 6, 8])
        ax.tick_params(labelsize=PT_TICK)
        ax.grid(axis="x", color="0.88", lw=0.4, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        print(f"  {T:3d} Ma: {n_tot} cells, CCD {c:.0f} m")

    axes[0].set_yticks(np.arange(1000, 6001, 1000))
    axes[0].set_yticklabels([f"{int(y/1000)}" for y in np.arange(1000, 6001, 1000)])
    axes[0].set_ylabel("Water depth (km)", fontsize=PT_LABEL)
    fig.supxlabel("Carbonate accumulation rate (m Myr$^{-1}$)", fontsize=PT_LABEL, y=0.035)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.92, bottom=0.145, wspace=0.14)

    bad = check_text_collisions(fig, list(axes))
    if bad:
        print("  ! text collisions:")
        for b in bad:
            print("     ", b)
    else:
        print("  no text collisions")

    for ext in ("png", "pdf"):
        out = FIGDIR / f"FigS5_accumulation_rate_vs_depth.{ext}"
        fig.savefig(out, dpi=400 if ext == "png" else None)
        print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
