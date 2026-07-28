#!/usr/bin/env python3
"""
Step 2 — Global CCD synthesis with an inter-basin dispersion envelope.

Combine the regional (Atlantic, Pacific, Indian) CCD curves into a single global
CCD curve over 0–52 Ma, and quantify the *inter-basin dispersion* envelope
(a physical heterogeneity metric, NOT an analytical error).

At each 1 Myr timestep the global CCD is the weighted mean of the basins that have
data, weighted by their ocean-area fractions (columns 7/8/9), renormalised over the
basins present. The dispersion is ``sigma(t) = sqrt(sum_i w_i (CCD_i - mean)^2)``,
defined only where >= 2 basins contribute.

This corrects the original ``2_CCD_averaging/1_CCD_averaging.py``, which used EQUAL
weights among available basins and ignored the fraction columns — a bug per the
authors. Set ``WEIGHTING = "equal"`` to reproduce that old behaviour.

Input:
    data/ccds_and_oceanbasin_fractions.xlsx  (header on row 2)
        col 0 Atlantic Age, 1 Atlantic CCD
        col 2 Pacific  Age, 3 Pacific  CCD
        col 4 Indian   Age, 5 Indian   CCD
        col 6 master age grid (1 Myr)
        col 7/8/9 Atlantic/Pacific/Indian AREA FRACTION  (see note below)

Output:
    outputs/step3_global_ccd/global_ccd_with_basin_dispersion_0-52Ma.txt
        Age_Ma, Global_CCD_m, CCD_minus_dispersion_m, CCD_plus_dispersion_m
        (integer values; commented header)

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


def synthesise(xlsx_path) -> pd.DataFrame:
    df = pd.read_excel(xlsx_path, header=1)
    for col in [1, 3, 5, 6, 7, 8, 9]:
        if col < df.shape[1]:
            df.iloc[:, col] = pd.to_numeric(df.iloc[:, col], errors="coerce")

    age = pd.to_numeric(df.iloc[:, 6], errors="coerce").to_numpy(float)
    ccd_cols = [df.iloc[:, 1].to_numpy(float),
                df.iloc[:, 3].to_numpy(float),
                df.iloc[:, 5].to_numpy(float)]
    frac_cols = None
    if WEIGHTING == "area" and df.shape[1] >= 10:
        frac_cols = [df.iloc[:, 7].to_numpy(float),
                     df.iloc[:, 8].to_numpy(float),
                     df.iloc[:, 9].to_numpy(float)]

    n = len(age)
    global_ccd = np.full(n, np.nan)
    dispersion = np.full(n, np.nan)

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
        if k >= 2:
            dispersion[i] = float(np.sqrt(np.sum(w * (cv - mu) ** 2)))

    mask = np.isfinite(age) & (age >= 0) & (age <= 52)
    age52, g52, d52 = age[mask], global_ccd[mask], dispersion[mask]

    lower = np.where(np.isfinite(d52), np.rint(g52 - d52), np.nan)
    upper = np.where(np.isfinite(d52), np.rint(g52 + d52), np.nan)

    return pd.DataFrame({
        "Age_Ma": np.rint(age52).astype(int),
        "Global_CCD_m": np.rint(g52).astype(int),
        "CCD_minus_dispersion_m": lower,
        "CCD_plus_dispersion_m": upper,
    })


def _plot(xlsx_path, out_df) -> None:
    """Diagnostic figure: the synthesised global CCD with its inter-basin
    dispersion envelope, over the three regional CCD curves."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ccdworkflow.io import save_matplotlib_figure

    from matplotlib.ticker import MultipleLocator
    df = pd.read_excel(xlsx_path, header=1)
    regionals = [("Atlantic", 0, 1, "#1f77b4"), ("Pacific", 2, 3, "#d62728"), ("Indian", 4, 5, "#2ca02c")]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for label, ca, cc, col in regionals:
        a = pd.to_numeric(df.iloc[:, ca], errors="coerce")
        c = pd.to_numeric(df.iloc[:, cc], errors="coerce")
        mask = a <= 52     # last 52 Ma only (match the published Fig. 1)
        ax.plot(a[mask], c[mask], color=col, lw=1.6, alpha=0.8, label=f"{label} CCD")
    age = out_df["Age_Ma"].to_numpy(float)
    if np.isfinite(out_df["CCD_minus_dispersion_m"]).any():
        ax.fill_between(age, out_df["CCD_minus_dispersion_m"], out_df["CCD_plus_dispersion_m"],
                        color="0.6", alpha=0.35, label="inter-basin dispersion")
    ax.plot(age, out_df["Global_CCD_m"], color="black", lw=2.6, label="Global CCD (area-weighted)")
    ax.set_xlim(52, 0)                        # 0–52 Ma only
    ax.set_xlabel("Age (Ma)", fontsize=13); ax.set_ylabel("CCD (m)", fontsize=13)
    ax.xaxis.set_major_locator(MultipleLocator(10))   # labelled ticks every 10 Myr
    ax.xaxis.set_minor_locator(MultipleLocator(5))    # frame tick marks every 5 Myr
    ax.tick_params(axis="x", which="both", direction="out")
    ax.tick_params(labelsize=11)
    ax.legend(frameon=False, fontsize=12)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step3_global_ccd_with_dispersion", dpi=300, step="step3")
    plt.close(fig)


def main() -> None:
    config.ensure_dirs()
    out_df = synthesise(config.CCD_FRACTIONS_XLSX)
    out = write_table(config.GLOBAL_CCD_0_52, out_df,
                      header=["Age_Ma", "Global_CCD_m",
                              "CCD_minus_dispersion_m", "CCD_plus_dispersion_m"],
                      float_format="%.0f")
    n_env = int(np.isfinite(out_df["CCD_minus_dispersion_m"]).sum())
    print(f"[step3] wrote {out}  ({len(out_df)} rows, dispersion defined for {n_env})")
    _plot(config.CCD_FRACTIONS_XLSX, out_df)


if __name__ == "__main__":
    main()
