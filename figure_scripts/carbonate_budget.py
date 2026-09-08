#!/usr/bin/env python3
"""
The carbonate budget: net carbonate volume gain per million years against the seafloor
area available above the CCD, at 1 Myr resolution from 170 Ma to the present.

Drawn as panel (d) of Figure 3 by make_fig3_maps.py, which calls draw_budget(). Run
this file directly to (re)write the data file on its own.

Left axis
    Net carbonate volume gain per Myr, V(t) - V(t+1), smoothed with a Gaussian of
    3 Myr full width at half maximum, where V(t) is the global
    compacted carbonate volume on the reconstructed seafloor at time t. This is a
    global integral, so unlike a cell-by-cell difference it is unaffected by plate
    motion carrying crust between grid cells. It is a NET quantity, deposition during
    the interval minus the carbonate carried down with whatever subducted in it, which
    is why the axis reads gain rather than added.

Right axis
    The area of seafloor lying above the CCD at each time, which is the area over
    which the model deposits carbonate at all.

Together they show the budget: the reservoir grows faster when more of the seafloor
sits above the compensation depth.

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

Output: Figures/Fig3_carbonate_budget.csv    (the plotted series)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

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
SMOOTH_FWHM_MYR = 3.0     # Gaussian full width at half maximum, in Myr
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


def smooth_gaussian(y: np.ndarray, fwhm_myr: float) -> np.ndarray:
    """Gaussian smoother of the given full width at half maximum, in Myr.

    The unsmoothed series carries single-Myr spikes that come from the reconstruction
    stepping crust in and out of the mask rather than from anything in the carbonate
    budget. At 1 Myr sampling a 3 Myr FWHM cuts the scatter between neighbouring points
    four-fold and the curvature ten-fold while keeping four fifths of the range. A 5 Myr
    width smooths harder but flattens real structure, notably the dip near 50 Ma and the
    late Neogene bump.
    Weights are renormalised over the samples actually available, so the ends are not
    pulled toward zero.
    """
    sigma = fwhm_myr / 2.3548
    half = max(1, int(round(3 * sigma)))
    k = np.arange(-half, half + 1)
    w = np.exp(-0.5 * (k / sigma) ** 2)
    out = np.empty_like(y, dtype=float)
    for i in range(len(y)):
        lo, hi = max(0, i - half), min(len(y), i + half + 1)
        ww = w[half - (i - lo): half + (hi - i)]
        out[i] = float(np.sum(y[lo:hi] * ww) / ww.sum())
    return out


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


def budget_series() -> dict:
    """The plotted series. Cheap enough (about 7 s) to recompute wherever it is drawn,
    so nothing has to trust a cached copy."""
    if not VOLCSV.exists():
        raise SystemExit(f"missing {VOLCSV}; run pipeline_carbon/02_carbonate_volume_stats.py")
    vol = pd.read_csv(VOLCSV).sort_values("Age_Ma")
    vol = vol[vol["Age_Ma"] <= TMAX]
    age = vol["Age_Ma"].to_numpy(float)
    V = vol["volume_1e6km3_mean"].to_numpy(float)
    # V(t) - V(t+1) is the net gain during the Myr ending at t. Ages ascend, so this is
    # -diff, reported at the midpoint of each interval.
    raw = -np.diff(V)
    gain_age = 0.5 * (age[:-1] + age[1:])
    gain = smooth_gaussian(raw, SMOOTH_FWHM_MYR)
    area = area_above_ccd([int(t) for t in age])
    ok = np.isfinite(area[:-1]) & np.isfinite(gain)
    r = float(np.corrcoef(area[:-1][ok], gain[ok])[0, 1])
    return dict(age=age, V=V, gain_age=gain_age, gain=gain, raw=raw, area=area, r=r)


def draw_budget(ax, s: dict, label_size: float = PT_LABEL, tick_size: float = PT_TICK):
    """Draw the budget onto an existing axes, returning the twinned right-hand axes.

    Shared so that the panel in Figure 3 of the paper and any standalone rendering come
    from one piece of code rather than two that can drift apart.
    """
    ax.plot(s["gain_age"], s["gain"], color=VOL_COLOUR, lw=1.6)
    ax.set_xlim(TMAX, 0)
    ax.set_xlabel("Age (Ma)", fontsize=label_size)
    ax.set_ylabel("Net carbonate volume gain\n(10$^6$ km$^3$ Myr$^{-1}$)",
                  fontsize=label_size, color=VOL_COLOUR)
    ax.tick_params(axis="y", labelcolor=VOL_COLOUR)
    ax.axhline(0, color="0.75", lw=0.6, zorder=0)

    bx = ax.twinx()
    bx.plot(s["age"], s["area"], color=AREA_COLOUR, lw=1.4, ls=(0, (5, 2)))
    bx.set_ylabel("Seafloor area above\nthe CCD (10$^6$ km$^2$)",
                  fontsize=label_size, color=AREA_COLOUR)
    bx.tick_params(axis="y", labelcolor=AREA_COLOUR)
    bx.set_ylim(0, None)
    for a in (ax, bx):
        a.tick_params(labelsize=tick_size)
        a.spines["top"].set_visible(False)
    ax.grid(axis="x", color="0.9", lw=0.4)
    ax.set_axisbelow(True)
    return bx


def write_csv(s: dict) -> Path:
    out_csv = FIGDIR / "Fig3_carbonate_budget.csv"
    pd.DataFrame({"Age_Ma": s["age"], "volume_1e6km3": s["V"],
                  "net_gain_1e6km3_per_myr_raw": np.concatenate([[np.nan], s["raw"]]),
                  "net_gain_1e6km3_per_myr_smoothed": np.concatenate([[np.nan], s["gain"]]),
                  "area_above_CCD_1e6km2": s["area"]}).to_csv(
        out_csv, index=False, float_format="%.4f")
    return out_csv


def text_collisions(fig, pairs) -> list:
    """Overlapping drawn text. `pairs` is a list of (axes, is_primary): x ticks are read
    only from the primary of a twinned pair, and ticks outside the view are skipped,
    since matplotlib keeps those and they would register as phantom overlaps."""
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()

    def live(axis, lim):
        lo, hi = min(lim), max(lim)
        return [l for loc, l in zip(axis.get_ticklocs(), axis.get_ticklabels())
                if lo <= loc <= hi and l.get_text().strip() and l.get_visible()]

    cand = []
    for a, primary in pairs:
        cand += list(a.texts) + [a.yaxis.label] + live(a.yaxis, a.get_ylim())
        if primary:
            cand += [a.xaxis.label] + live(a.xaxis, a.get_xlim())
    items = []
    for t in cand:
        if not (t.get_text().strip() and t.get_visible()):
            continue
        e = t.get_window_extent(renderer=rr)
        if e.width > 0 and e.height > 0:
            items.append((t, e))
    bad = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a_, b_ = items[i][1], items[j][1]
            if (a_.overlaps(b_) and min(a_.x1, b_.x1) - max(a_.x0, b_.x0) > 1.0
                    and min(a_.y1, b_.y1) - max(a_.y0, b_.y0) > 1.0):
                bad.append(f"{items[i][0].get_text()!r} / {items[j][0].get_text()!r}")
    return bad


def main() -> None:
    """Standalone use writes only the data file. The figure itself is now panel (d) of
    Figure 3 in the paper, drawn by make_fig3_maps.py through draw_budget()."""
    s = budget_series()
    print(f"  correlation of net carbonate gain with area above the CCD: r = {s['r']:.2f}")
    print(f"wrote {write_csv(s)}")


if __name__ == "__main__":
    main()
