#!/usr/bin/env python3
"""
Figure S5 - carbonate added per million years, against the seafloor area available
above the CCD, at 1 Myr resolution from 170 Ma to the present.

Left axis
    Net carbonate volume added per Myr, V(t) - V(t+1), where V(t) is the global
    compacted carbonate volume on the reconstructed seafloor at time t. This is a
    global integral, so unlike a cell-by-cell difference it is unaffected by plate
    motion carrying crust between grid cells. It is a NET quantity: deposition during
    the interval minus the carbonate carried down with whatever subducted in it.

Right axis
    The area of seafloor lying above the CCD at each time, which is the area over
    which the model deposits carbonate at all.

Together they show the budget: carbonate accumulates faster when more of the
seafloor sits above the compensation depth.

Definition of the seafloor area
-------------------------------
Cells where the carbonate-thickness grid and the paleobathymetry are both defined
and the bathymetry is below sea level, that is, the model's oceanic crust. Areas are
true spherical cell areas, not cell counts. "Above the CCD" means the reconstructed
bathymetry at that time is shallower than the CCD from this study at that time. No
depth threshold is applied, so this is ocean floor rather than abyssal plain
specifically; the total is 307 x 10^6 km^2 at the present day.

Inputs:
    steps/step9_carbonate_volume_analysis/output/dm2026_volume_stats/
        carbonate_volume_area_through_time.csv        (global volume per 1 Myr)
    steps/step8_carbonate_sediment_thickness/carbonate_sed_thickness_DM2026/
        compacted_sediment_thickness_0.25_{t}.nc      (ocean mask)
    steps/step8_carbonate_sediment_thickness/input_grids/
        Alfonso2024_pybacktrack_merged_paleobathymetry/paleobathymetry_{t}Ma.nc
    Figures/CCD_hybrid_DM2026.txt

Output: Figures/FigS5_carbonate_budget.{png,pdf}
        Figures/FigS5_carbonate_budget.csv    (the plotted series)
"""
from __future__ import annotations

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
CARB = STEP8 / "carbonate_sed_thickness_DM2026"
BATHY = STEP8 / "input_grids/Alfonso2024_pybacktrack_merged_paleobathymetry"
VOLCSV = (CW / "steps/step9_carbonate_volume_analysis/output/dm2026_volume_stats"
             / "carbonate_volume_area_through_time.csv")
HYBRID = FIGDIR / "CCD_hybrid_DM2026.txt"

TMAX = 170
R_EARTH_M = 6371000.0
VOL_COLOUR, AREA_COLOUR = "#1f6f8b", "#b3202c"
PT_TICK, PT_LABEL, PT_ANNOT = 7.5, 8.5, 7.5


def _pm180(da: xr.DataArray) -> xr.DataArray:
    if float(da["lon"].max()) > 180.0:
        da = da.assign_coords(lon=((da.lon + 180.0) % 360.0) - 180.0).sortby("lon")
        _, keep = np.unique(da["lon"].values, return_index=True)
        da = da.isel(lon=np.sort(keep))
    return da


def _cell_area_m2(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """True spherical area of each cell in a latitude row."""
    dlat = abs(lat[1] - lat[0])
    dlon = abs(lon[1] - lon[0])
    edges = np.deg2rad(np.clip(np.concatenate(
        [[lat[0] - dlat / 2], (lat[:-1] + lat[1:]) / 2, [lat[-1] + dlat / 2]]), -90, 90))
    return R_EARTH_M ** 2 * np.deg2rad(dlon) * np.abs(np.sin(edges[1:]) - np.sin(edges[:-1]))


def ccd_curve():
    h = pd.read_csv(HYBRID, sep=r"\s+", header=None, names=["age", "ccd", "lo", "hi"])
    return lambda t: float(np.interp(t, h["age"], h["ccd"]))     # negative down


def area_above_ccd(times) -> np.ndarray:
    """Seafloor area shallower than the CCD, in 1e6 km^2, at each time."""
    ccd_of = ccd_curve()
    out = np.full(len(times), np.nan)
    band = None
    for i, T in enumerate(times):
        cp = CARB / f"compacted_sediment_thickness_0.25_{T}.nc"
        bp = BATHY / f"paleobathymetry_{T}Ma.nc"
        if not (cp.exists() and bp.exists()):
            continue
        carb = xr.open_dataset(cp)["z"]
        if band is None:
            band = _cell_area_m2(carb.lat.values, carb.lon.values)
            A = np.broadcast_to(band[:, None], carb.shape)
        B = _pm180(xr.open_dataset(bp)["z"]).interp(
            lon=carb.lon, lat=carb.lat, method="linear").values
        ocean = np.isfinite(B) & np.isfinite(carb.values) & (B < 0.0)
        out[i] = float(A[ocean & (B > ccd_of(T))].sum()) / 1e12      # m^2 -> 1e6 km^2
        if T % 20 == 0:
            print(f"  {T:3d} Ma: {out[i]:6.1f} x 1e6 km2 above the CCD")
    return out


def main() -> None:
    if not VOLCSV.exists():
        raise SystemExit(f"missing {VOLCSV}; run pipeline_carbon/02_carbonate_volume_stats.py")
    vol = pd.read_csv(VOLCSV).sort_values("Age_Ma")
    vol = vol[vol["Age_Ma"] <= TMAX]
    age = vol["Age_Ma"].to_numpy(float)
    V = vol["volume_1e6km3_mean"].to_numpy(float)
    # V(t) - V(t+1): carbonate added during the Myr ending at t. Ages ascend, so this
    # is -diff. Reported at the midpoint of each interval.
    added = -np.diff(V)
    added_age = 0.5 * (age[:-1] + age[1:])

    times = [int(t) for t in age]
    area = area_above_ccd(times)

    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    ax.plot(added_age, added, color=VOL_COLOUR, lw=1.4)
    ax.set_xlim(TMAX, 0)
    ax.set_xlabel("Age (Ma)", fontsize=PT_LABEL)
    ax.set_ylabel("Carbonate added (10$^6$ km$^3$ Myr$^{-1}$)", fontsize=PT_LABEL,
                  color=VOL_COLOUR)
    ax.tick_params(axis="y", labelcolor=VOL_COLOUR)
    ax.axhline(0, color="0.75", lw=0.6, zorder=0)

    bx = ax.twinx()
    bx.plot(age, area, color=AREA_COLOUR, lw=1.4, ls=(0, (5, 2)))
    bx.set_ylabel("Seafloor above the CCD (10$^6$ km$^2$)", fontsize=PT_LABEL,
                  color=AREA_COLOUR)
    bx.tick_params(axis="y", labelcolor=AREA_COLOUR)
    bx.set_ylim(0, None)

    for a in (ax, bx):
        a.tick_params(labelsize=PT_TICK)
    ax.grid(axis="x", color="0.9", lw=0.4)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    bx.spines["top"].set_visible(False)

    ok = np.isfinite(area[:-1]) & np.isfinite(added)
    r = float(np.corrcoef(area[:-1][ok], added[ok])[0, 1])
    print(f"  correlation of carbonate added with area above the CCD: r = {r:.2f}")

    out_csv = FIGDIR / "FigS5_carbonate_budget.csv"
    pd.DataFrame({"Age_Ma": age, "volume_1e6km3": V,
                  "area_above_CCD_1e6km2": area}).to_csv(out_csv, index=False,
                                                         float_format="%.4f")
    print(f"wrote {out_csv}")

    fig.subplots_adjust(left=0.105, right=0.885, top=0.965, bottom=0.165)
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    def live_ticklabels(axis, lo, hi):
        """Tick labels that are actually drawn: matplotlib keeps the ones for ticks
        outside the view limits, and those otherwise register as phantom collisions."""
        lo, hi = min(lo, hi), max(lo, hi)
        return [lab for loc, lab in zip(axis.get_ticklocs(), axis.get_ticklabels())
                if lo <= loc <= hi and lab.get_text().strip() and lab.get_visible()]

    # On a twin axis both axes carry the same x ticks in the same place, so take the
    # x ticks from the primary axis only; otherwise every one collides with its twin.
    cand = list(ax.texts) + [ax.xaxis.label, ax.yaxis.label, bx.yaxis.label] + list(bx.texts)
    cand += live_ticklabels(ax.xaxis, *ax.get_xlim())
    cand += live_ticklabels(ax.yaxis, *ax.get_ylim())
    cand += live_ticklabels(bx.yaxis, *bx.get_ylim())
    items = []
    for t in cand:
        if not (t.get_text().strip() and t.get_visible()):
            continue
        e = t.get_window_extent(renderer=rr)
        if e.width > 0 and e.height > 0:
            items.append((t, e))
    bad = [f"{items[i][0].get_text()!r} / {items[j][0].get_text()!r}"
           for i in range(len(items)) for j in range(i + 1, len(items))
           if items[i][1].overlaps(items[j][1])
           and min(items[i][1].x1, items[j][1].x1) - max(items[i][1].x0, items[j][1].x0) > 1.0
           and min(items[i][1].y1, items[j][1].y1) - max(items[i][1].y0, items[j][1].y0) > 1.0]
    print("  text collisions: " + (", ".join(bad) if bad else "none"))

    for ext in ("png", "pdf"):
        out = FIGDIR / f"FigS5_carbonate_budget.{ext}"
        fig.savefig(out, dpi=400 if ext == "png" else None)
        print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
