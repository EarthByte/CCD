#!/usr/bin/env python3
"""
Figure S5 - the model's carbonate accumulation rate against water depth, and the
share of abyssal seafloor on which carbonate accumulates, from 80 Ma to the present.

What is plotted
---------------
(a) The instantaneous decompacted carbonate accumulation rate as a function of water
    depth at 80, 60, 40, 20 and 0 Ma. This is the rate the model applies at each
    reconstruction time,

        rate = max_rate(t) * min(bathymetry - CCD, 300 m) / 300 m * latitude_factor

    evaluated on the paleobathymetry grid for that time, with the CCD from this study.
    Carbonate accumulates at the full rate more than 300 m above the CCD, tapers to
    zero through the dissolution zone, and is zero below it. The curve is the median
    over all abyssal cells in each 100 m depth interval between 2 and 6 km; shallower
    than that there is too little abyssal seafloor for a stable median. At a given
    depth the only source of spread is the latitude factor, which ramps in between
    34 and 23 Ma.

(b) The fraction of abyssal seafloor lying above the CCD, area-weighted, every 5 Myr.
    Unlike the profiles in (a), this is not prescribed: it follows from the
    reconstructed bathymetry meeting the reconstructed CCD, and it is the quantity
    that actually changes through time.

Why not a rate from differencing thickness grids
------------------------------------------------
Thickness grids at adjacent times cannot be differenced. The grids are reconstructed
into palaeo-coordinates, and at 0.25 degrees a plate moving 20-100 km/Myr crosses one
to four cells per Myr, so the difference is dominated by lateral advection of the
thickness gradient rather than by deposition: 36-40% of cells return a negative rate,
which the model cannot produce, and after clipping those to zero roughly a tenth of
the remainder still exceed the model's own 18 m/Myr ceiling.

Why not violins
---------------
The instantaneous rate is a deterministic function of depth below the CCD, so it has
almost no distribution to show: before 34 Ma every cell more than 300 m above the CCD
sits exactly at the ceiling, and after the latitude factor switches on 64-76% still
do. A violin of that is a spike, so the profile is drawn as a median line with a
percentile band instead.

Inputs:
    steps/step8_carbonate_sediment_thickness/input_grids/
        Alfonso2024_pybacktrack_merged_paleobathymetry/paleobathymetry_{T}Ma.nc
    steps/step8_carbonate_sediment_thickness/input_data/sed_rate_best.txt
    steps/step8_carbonate_sediment_thickness/latitude_factor.py
    steps/step8_carbonate_sediment_thickness/carbonate_sed_thickness_DM2026/  (ocean mask)
    Figures/CCD_hybrid_DM2026.txt

Output: Figures/FigS5_accumulation_rate_vs_depth.{png,pdf}
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
BATHY = STEP8 / "input_grids/Alfonso2024_pybacktrack_merged_paleobathymetry"
CARB = STEP8 / "carbonate_sed_thickness_DM2026"
HYBRID = FIGDIR / "CCD_hybrid_DM2026.txt"

TIMES = (80, 60, 40, 20, 0)
AREA_TIMES = range(0, 81, 5)
DISSOLUTION_M = 300.0        # CCD_DISSOLUTION_DISTANCE in carbonate_sediment_thickness.py
DEPTH_STEP = 100.0
DEPTH_MIN = 2000.0        # shallower than this there is too little abyssal
DEPTH_MAX = 6000.0        # seafloor for the median to be stable
MIN_CELLS = 1000
PT_TICK, PT_LABEL, PT_PANEL, PT_ANNOT = 7.5, 8.5, 10.0, 7.5


def _latitude_factor():
    spec = importlib.util.spec_from_file_location("lf", STEP8 / "latitude_factor.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["lf"] = m
    spec.loader.exec_module(m)
    return m.latitude_factor


def _max_rate_curve():
    """sed_rate_best.txt is the decompacted maximum rate in cm/kyr; 10x gives m/Myr."""
    rows = []
    for line in open(STEP8 / "input_data/sed_rate_best.txt"):
        s = line.split()
        if len(s) >= 2:
            try:
                rows.append((float(s[0]), float(s[1])))
            except ValueError:
                pass
    a = np.array(sorted(rows))
    return lambda t: float(np.interp(t, a[:, 0], a[:, 1])) * 10.0


def _ccd_curve():
    h = pd.read_csv(HYBRID, sep=r"\s+", header=None, names=["age", "ccd", "lo", "hi"])
    return lambda t: float(np.interp(t, h["age"], h["ccd"]))     # negative down


def _pm180(da: xr.DataArray) -> xr.DataArray:
    if float(da["lon"].max()) > 180.0:
        da = da.assign_coords(lon=((da.lon + 180.0) % 360.0) - 180.0).sortby("lon")
        _, keep = np.unique(da["lon"].values, return_index=True)
        da = da.isel(lon=np.sort(keep))
    return da


def fields(T: int, lat_factor, max_rate, ccd_of):
    """Depth (m, positive down), rate (m/Myr), and cell area weight, for abyssal cells."""
    carb = xr.open_dataset(CARB / f"compacted_sediment_thickness_0.25_{T}.nc")["z"]
    bath = _pm180(xr.open_dataset(BATHY / f"paleobathymetry_{T}Ma.nc")["z"])
    B = bath.interp(lon=carb.lon, lat=carb.lat, method="linear").values
    lat = carb.lat.values
    LAT = np.broadcast_to(lat[:, None], B.shape)
    ok = np.isfinite(B) & np.isfinite(carb.values) & (B < 0.0)
    uniq = np.unique(lat)
    lut = {y: lat_factor(float(y), T) for y in uniq}
    LF = np.vectorize(lut.get)(LAT)
    ccd = ccd_of(T)
    rate = np.where(B > ccd, max_rate(T) * np.clip(B - ccd, 0.0, DISSOLUTION_M)
                    / DISSOLUTION_M * LF, 0.0)
    w = np.broadcast_to(np.cos(np.deg2rad(lat))[:, None], B.shape)
    return -B[ok], rate[ok], w[ok], ccd


def main() -> None:
    lat_factor, max_rate, ccd_of = _latitude_factor(), _max_rate_curve(), _ccd_curve()
    edges = np.arange(DEPTH_MIN, DEPTH_MAX + DEPTH_STEP, DEPTH_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])
    colours = plt.get_cmap("viridis")(np.linspace(0.08, 0.82, len(TIMES)))

    fig = plt.figure(figsize=(7.48, 3.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.30,
                          left=0.085, right=0.985, top=0.95, bottom=0.155)
    ax, bx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    for col, T in zip(colours, TIMES):
        depth, rate, _, ccd = fields(T, lat_factor, max_rate, ccd_of)
        med = np.full(len(centres), np.nan)
        for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
            k = (depth >= a) & (depth < b)
            if k.sum() >= MIN_CELLS:
                med[i] = np.median(rate[k])
        ax.plot(med, centres, color=col, lw=1.6, label=f"{T} Ma")
        ax.plot([0.0], [-ccd], marker="v", ms=4.5, color=col, clip_on=False)
        print(f"  {T:3d} Ma: CCD {-ccd:.0f} m, max rate {max_rate(T):.1f} m/Myr")

    ax.set_xlim(0, 19.5)
    ax.set_ylim(DEPTH_MAX, DEPTH_MIN)
    ax.set_yticks(np.arange(DEPTH_MIN, 6001, 500))
    ax.set_yticklabels([f"{y/1000:.1f}" for y in np.arange(DEPTH_MIN, 6001, 500)])
    ax.set_xlabel("Carbonate accumulation rate (m Myr$^{-1}$)", fontsize=PT_LABEL)
    ax.set_ylabel("Water depth (km)", fontsize=PT_LABEL)
    ax.legend(frameon=False, fontsize=PT_ANNOT, loc="lower right", handlelength=1.4,
              title="triangles mark the CCD", title_fontsize=PT_ANNOT)
    ax.get_legend().get_title().set_color("0.35")

    frac = []
    for T in AREA_TIMES:
        depth, rate, w, ccd = fields(T, lat_factor, max_rate, ccd_of)
        frac.append(100.0 * w[depth < -ccd].sum() / w.sum())
    bx.plot(list(AREA_TIMES), frac, color="0.2", lw=1.6)
    for col, T in zip(colours, TIMES):
        bx.plot([T], [frac[list(AREA_TIMES).index(T)]], marker="o", ms=5, color=col, zorder=5)
    bx.set_xlim(80, 0)
    bx.set_ylim(0, 100)
    bx.set_xlabel("Age (Ma)", fontsize=PT_LABEL)
    bx.set_ylabel("Abyssal seafloor above the CCD (%)", fontsize=PT_LABEL)
    print("  above-CCD area fraction (%): " +
          ", ".join(f"{t}:{f:.0f}" for t, f in zip(AREA_TIMES, frac)))

    for a, letter in ((ax, "a"), (bx, "b")):
        a.tick_params(labelsize=PT_TICK)
        a.grid(color="0.9", lw=0.4)
        a.set_axisbelow(True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        a.text(-0.02, 1.02, letter, transform=a.transAxes, fontsize=PT_PANEL,
               fontweight="bold", ha="left", va="bottom")

    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = [(t, t.get_window_extent(renderer=r)) for a in (ax, bx)
             for t in a.texts + a.get_xticklabels() + a.get_yticklabels()
             + [a.xaxis.label, a.yaxis.label] if t.get_text().strip()]
    bad = [f"{items[i][0].get_text()!r} / {items[j][0].get_text()!r}"
           for i in range(len(items)) for j in range(i + 1, len(items))
           if items[i][1].overlaps(items[j][1])
           and min(items[i][1].x1, items[j][1].x1) - max(items[i][1].x0, items[j][1].x0) > 1.0
           and min(items[i][1].y1, items[j][1].y1) - max(items[i][1].y0, items[j][1].y0) > 1.0]
    print("  text collisions: " + (", ".join(bad) if bad else "none"))

    for ext in ("png", "pdf"):
        out = FIGDIR / f"FigS5_accumulation_rate_vs_depth.{ext}"
        fig.savefig(out, dpi=400 if ext == "png" else None)
        print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
