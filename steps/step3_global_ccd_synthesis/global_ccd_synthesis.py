#!/usr/bin/env python3
"""
Step 3 — Global CCD synthesis with an inter-basin dispersion envelope.

Combine the regional (Atlantic, Pacific, Indian) CCD curves into a single global
CCD curve over 0–52 Ma, and quantify the *inter-basin dispersion* envelope
(a physical heterogeneity metric, NOT an analytical error).

At each 1 Myr timestep the global CCD is the weighted mean of the basins that have
data, weighted by their ocean-area fractions (columns 7/8/9), renormalised over the
basins present. The plotted envelope is the FULL INTER-BASIN RANGE, i.e. the
minimum and maximum of the contributing basin CCDs, defined only where >= 2
basins contribute.

Why not the weighted standard deviation
---------------------------------------
The envelope used to be ``sigma(t) = sqrt(sum_i w_i (CCD_i - mean)^2)``. With uneven
area weights that band is NOT bounded by the basin curves it is drawn over: for two
basins with weights p and 1-p, sigma = sqrt(p(1-p))*|delta| while the distance from
the weighted mean to the nearer curve is only min(p,1-p)*|delta|, and
sqrt(p(1-p)) > min(p,1-p) for every p != 0.5. In this dataset the +-1 sigma band fell
outside the basin range at 51 of 53 timesteps, by up to 292 m, which is what made
Fig. 1 look wrong. The weighted sigma is still written to
``outputs/inter_basin_stats_0-52Ma.txt`` for reference.

This corrects the original ``2_CCD_averaging/1_CCD_averaging.py``, which used EQUAL
weights among available basins and ignored the fraction columns — a bug per the
authors. Set ``WEIGHTING = "equal"`` to reproduce that old behaviour.

Inputs:
    steps/step1_timescale_conversion/outputs/{ATL,PAC,IND}_CCD_GTS2020_1my.txt
        the three regional CCD curves, all normalised to GTS2020 and resampled
        onto a 1 Myr grid by step 0
    data/ccds_and_oceanbasin_fractions.xlsx  (header on row 2)
        col 6 master age grid (1 Myr)
        col 7/8/9 Atlantic/Pacific/Indian AREA FRACTION  (see note below)
        The regional CCD columns 0-5 of this spreadsheet are NO LONGER read here:
        they carry the published curves on their original, mutually inconsistent
        timescales. Step 0 converts them; this step consumes the converted files.

Output:
    outputs/step3_global_ccd/global_ccd_with_basin_dispersion_0-52Ma.txt
        Age_Ma, Global_CCD_m, CCD_basin_min_m, CCD_basin_max_m
        (integer values; commented header)
    outputs/inter_basin_stats_0-52Ma.txt
        Age_Ma, N_basins, Global_CCD_m, Basin_min_m, Basin_max_m, Weighted_SD_m

Weighting note
--------------
The paper text describes an *area-weighted* combination, and the xlsx now carries
basin-area fraction columns (7/8/9). The canonical script that produced the
published curve nonetheless used EQUAL weights among available basins. That
behaviour is reproduced here (``WEIGHTING = "equal"``) so numbers match. Set
``WEIGHTING = "area"`` to instead weight by the fraction columns. See README
"Known discrepancies".
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ccdworkflow import config
from ccdworkflow.io import write_table

WEIGHTING = "area"    # "area" (area-weighted, intended) or "equal" (old buggy behaviour)


BASINS = ("Atlantic", "Pacific", "Indian")


def load_regional_ccds() -> dict:
    """The GTS2020-normalised regional CCD curves written by step 0."""
    out = {}
    for name in BASINS:
        path = config.GTS2020_REGIONAL_CCD[name]
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing - run step 0 (timescale conversion) first."
            )
        d = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                        names=["Age_Ma", "CCD_m"])
        a = pd.to_numeric(d["Age_Ma"], errors="coerce").to_numpy(float)
        v = pd.to_numeric(d["CCD_m"], errors="coerce").to_numpy(float)
        m = np.isfinite(a) & np.isfinite(v)
        out[name] = (a[m], v[m])
    return out


def _onto(grid, ages, values):
    """Sample a regional curve on the master grid, NaN outside its own age span."""
    y = np.interp(grid, ages, values, left=np.nan, right=np.nan)
    return np.where((grid >= ages.min()) & (grid <= ages.max()), y, np.nan)


def synthesise(xlsx_path) -> pd.DataFrame:
    df = pd.read_excel(xlsx_path, header=1)
    for col in [6, 7, 8, 9]:
        if col < df.shape[1]:
            df.iloc[:, col] = pd.to_numeric(df.iloc[:, col], errors="coerce")

    age = pd.to_numeric(df.iloc[:, 6], errors="coerce").to_numpy(float)
    reg = load_regional_ccds()
    ccd_cols = [_onto(age, *reg[name]) for name in BASINS]
    frac_cols = None
    if WEIGHTING == "area" and df.shape[1] >= 10:
        frac_cols = [df.iloc[:, 7].to_numpy(float),
                     df.iloc[:, 8].to_numpy(float),
                     df.iloc[:, 9].to_numpy(float)]

    n = len(age)
    global_ccd = np.full(n, np.nan)
    band_lo = np.full(n, np.nan)      # minimum of the contributing basins
    band_hi = np.full(n, np.nan)      # maximum of the contributing basins
    wsd = np.full(n, np.nan)          # weighted SD, kept for reference only
    nbas = np.zeros(n, dtype=int)

    for i in range(n):
        ccds = np.array([c[i] for c in ccd_cols], dtype=float)
        valid = np.isfinite(ccds)
        k = int(valid.sum())
        if k == 0:
            continue
        cv = ccds[valid]
        if frac_cols is not None:
            w = np.array([frac_cols[j][i] for j in range(3)], dtype=float)[valid]
            w = np.nan_to_num(w, nan=0.0)
            w = w / w.sum() if w.sum() > 0 else np.ones(k) / k
        else:
            w = np.ones(k) / k
        mu = float(np.sum(w * cv))
        global_ccd[i] = mu
        nbas[i] = k
        if k >= 2:
            band_lo[i] = float(cv.min())
            band_hi[i] = float(cv.max())
            wsd[i] = float(np.sqrt(np.sum(w * (cv - mu) ** 2)))

    mask = np.isfinite(age) & (age >= 0) & (age <= 52)
    age52, g52 = age[mask], global_ccd[mask]
    lo52, hi52, sd52, k52 = band_lo[mask], band_hi[mask], wsd[mask], nbas[mask]

    main = pd.DataFrame({
        "Age_Ma": np.rint(age52).astype(int),
        "Global_CCD_m": np.rint(g52).astype(int),
        "CCD_basin_min_m": np.where(np.isfinite(lo52), np.rint(lo52), np.nan),
        "CCD_basin_max_m": np.where(np.isfinite(hi52), np.rint(hi52), np.nan),
    })
    stats = pd.DataFrame({
        "Age_Ma": np.rint(age52).astype(int),
        "N_basins": k52,
        "Global_CCD_m": np.rint(g52).astype(int),
        "Basin_min_m": np.where(np.isfinite(lo52), np.rint(lo52), np.nan),
        "Basin_max_m": np.where(np.isfinite(hi52), np.rint(hi52), np.nan),
        "Weighted_SD_m": np.where(np.isfinite(sd52), np.rint(sd52), np.nan),
    })
    return main, stats


def _plot(out_df) -> None:
    """Diagnostic figure: the synthesised global CCD with its inter-basin
    dispersion envelope, over the three regional CCD curves."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ccdworkflow.io import save_matplotlib_figure

    from matplotlib.ticker import MultipleLocator
    reg = load_regional_ccds()
    colours = {"Atlantic": "#1f77b4", "Pacific": "#d62728", "Indian": "#2ca02c"}
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for label in BASINS:
        a, c = reg[label]
        mask = a <= 52     # last 52 Ma only (match the published Fig. 1)
        ax.plot(a[mask], c[mask], color=colours[label], lw=1.6, alpha=0.8,
                label=f"{label} CCD")
    age = out_df["Age_Ma"].to_numpy(float)
    if np.isfinite(out_df["CCD_basin_min_m"]).any():
        ax.fill_between(age, out_df["CCD_basin_min_m"], out_df["CCD_basin_max_m"],
                        color="0.6", alpha=0.35, label="inter-basin range")
    ax.plot(age, out_df["Global_CCD_m"], color="black", lw=2.6, label="Global CCD (area-weighted)")
    ax.set_xlim(52, 0)                        # 0–52 Ma only
    ax.set_xlabel("Age (Ma)", fontsize=13); ax.set_ylabel("CCD (m)", fontsize=13)
    ax.xaxis.set_major_locator(MultipleLocator(10))   # labelled ticks every 10 Myr
    ax.xaxis.set_minor_locator(MultipleLocator(5))    # frame tick marks every 5 Myr
    ax.tick_params(axis="x", which="both", direction="out")
    ax.tick_params(labelsize=11)
    ax.legend(frameon=False, fontsize=12)
    fig.tight_layout()
    # This diagnostic IS Figure 1 of the paper, so it is written straight into the
    # paper's figure folder as well - no manual copy step.
    save_matplotlib_figure(fig, "step3_global_ccd_with_dispersion", dpi=300, step="step3",
                           paper_name="Fig1_regional_global_ccd")
    plt.close(fig)


def main() -> None:
    config.ensure_dirs()
    out_df, stats_df = synthesise(config.CCD_FRACTIONS_XLSX)
    out = write_table(config.GLOBAL_CCD_0_52, out_df,
                      header=["Age_Ma", "Global_CCD_m",
                              "CCD_basin_min_m", "CCD_basin_max_m"],
                      float_format="%.0f")
    write_table(config.STEP3_DIR / "inter_basin_stats_0-52Ma.txt", stats_df,
                header=list(stats_df.columns), float_format="%.0f")
    n_env = int(np.isfinite(out_df["CCD_basin_min_m"]).sum())
    print(f"[step3] wrote {out}  ({len(out_df)} rows, dispersion defined for {n_env})")
    _plot(out_df)


if __name__ == "__main__":
    main()
