#!/usr/bin/env python3
"""
Step 8 (batch, DM2026-only) — carbonate volume & area through time for the
MIN / MEAN / MAX carbonate-thickness grids.

The original step-8 scripts (02–06) compute DM2026-vs-BW1991 *comparison*
statistics; with the BW1991 comparison dropped, the natural DM2026-only statistic
is the global preserved carbonate volume and carbonate-covered sea-floor area
through time, with the min–max grids giving an uncertainty envelope.

For each scenario and each 1-Myr grid it computes:
    volume  = sum(thickness_m * cell_area_m2)          [m^3]  -> reported in 1e6 km^3
    area    = sum(cell_area_m2 where thickness finite) [m^2]  -> reported in 1e6 km^2

Outputs (config STEP8 output dir):
    carbonate_volume_area_through_time.csv
    carbonate_volume_through_time.{png,pdf}   (mean with min–max envelope)

Reads compacted-thickness NetCDFs written by step 7
(``compacted_sediment_thickness_0.25_<t>.nc``, variable ``z``). Runs anywhere the
grids exist and xarray/netCDF4 are installed.
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

import os
HERE = Path(__file__).resolve().parent
# steps/ sits beside pipeline_carbon/ inside the workflow, not one level further up.
# The runner exports STEPS_ROOT so the old default never bit under run_all.sh, but a
# standalone run silently wrote its outputs to a stray folder outside the workflow.
_STEPS_ROOT = Path(os.environ.get("STEPS_ROOT", HERE.parent / "steps"))
if not _STEPS_ROOT.is_dir():
    raise SystemExit(f"steps root not found: {_STEPS_ROOT} (set STEPS_ROOT to override)")
STEP7 = _STEPS_ROOT / "step8_carbonate_sediment_thickness"
OUT_DIR = _STEPS_ROOT / "step9_carbonate_volume_analysis" / "output" / "dm2026_volume_stats"

SCENARIOS = {
    "mean": STEP7 / "carbonate_sed_thickness_DM2026",
    "min":  STEP7 / "carbonate_sed_thickness_min_DM2026",
    "max":  STEP7 / "carbonate_sed_thickness_max_DM2026",
}
FMT = "compacted_sediment_thickness_0.25_{t}.nc"
R = 6_371_000.0   # m


def cell_area_m2(lat, lon):
    lat = np.asarray(lat, float); lon = np.asarray(lon, float)
    dlon = np.deg2rad(np.abs(np.median(np.diff(lon))))
    dlat = np.abs(np.median(np.diff(lat)))
    lat_edges = np.deg2rad(lat)
    band = R * R * dlon * np.abs(np.cos(lat_edges)) * np.deg2rad(dlat)   # per-cell along a lat row
    return band  # length nlat; broadcast across lon


def stats_for_grid(path: Path):
    ds = xr.open_dataset(path)
    z = ds["z"]
    lat = ds["lat"].values; lon = ds["lon"].values
    band = cell_area_m2(lat, lon)[:, None]          # (nlat,1) m^2 per cell
    area = np.broadcast_to(band, z.shape)
    zv = z.values
    finite = np.isfinite(zv)
    vol = float(np.nansum(np.where(finite, zv, 0.0) * area))       # m^3
    cov = float(np.nansum(np.where(finite, area, 0.0)))            # m^2
    ds.close()
    return vol, cov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmin", type=int, default=0)
    ap.add_argument("--tmax", type=int, default=170)
    ap.add_argument("--outdir", default=str(OUT_DIR))
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for t in range(args.tmin, args.tmax + 1):
        row = {"Age_Ma": t}
        for name, d in SCENARIOS.items():
            p = d / FMT.format(t=t)
            if p.exists():
                vol, cov = stats_for_grid(p)
                # 1e6 km^3 = 1e15 m^3 (1 km^3 = 1e9 m^3); 1e6 km^2 = 1e12 m^2 (1 km^2 = 1e6 m^2).
                # Both were previously divided by 1e18, understating volume 1000-fold and
                # area a million-fold. The curve shapes were unaffected, the axis numbers
                # were not.
                row[f"volume_1e6km3_{name}"] = vol / 1e15      # m^3 -> 1e6 km^3
                row[f"area_1e6km2_{name}"] = cov / 1e12        # m^2 -> 1e6 km^2
            else:
                row[f"volume_1e6km3_{name}"] = np.nan
                row[f"area_1e6km2_{name}"] = np.nan
        rows.append(row)
    df = pd.DataFrame(rows)
    csv = outdir / "carbonate_volume_area_through_time.csv"
    df.to_csv(csv, index=False, float_format="%.4f")
    print(f"[stats] wrote {csv}  ({len(df)} times)")

    if df["volume_1e6km3_mean"].notna().any():
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.fill_between(df["Age_Ma"], df["volume_1e6km3_min"], df["volume_1e6km3_max"],
                        color="lightsteelblue", alpha=0.7, label="min–max CCD envelope")
        ax.plot(df["Age_Ma"], df["volume_1e6km3_mean"], color="navy", lw=2, label="mean CCD")
        ax.invert_xaxis()
        ax.set_xlabel("Age (Ma)"); ax.set_ylabel("Preserved deep-sea carbonate volume (10$^6$ km$^3$)")
        ax.set_title("Global preserved deep-sea carbonate volume through time (DM2026)")
        ax.legend(frameon=False)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(outdir / f"carbonate_volume_through_time.{ext}", dpi=200)
        plt.close(fig)
        print(f"[stats] wrote carbonate_volume_through_time.png/pdf")


if __name__ == "__main__":
    main()
