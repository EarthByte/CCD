#!/usr/bin/env python3
"""
Sea-level vs CCD time-series regression diagnostics (0–52 Ma).

Faithful, lightly-refactored port of the original
``5_regression/ccd_sealevel_timeseries_regression.py``. Reports OLS, reduced major
axis (RMA, preferred), orthogonal distance regression (ODR) and Theil–Sen robust
fits; raw / detrended / first-difference tests; AR(1)-adjusted p-values; a
moving-block bootstrap; and a lag test. CCD is treated as positive-down depth
(``CCD_m = -Global_CCD_m``), so a negative slope means higher sea level ->
shallower CCD.

Call :func:`run` from step 5, or run standalone.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import odr, stats


@dataclass
class FitResult:
    method: str
    slope: float
    intercept: float
    r: float
    r2: float
    p_naive: float
    n: int
    neff_ar1: float
    p_ar1: float
    notes: str = ""


def _read2(path: Path, names):
    df = pd.read_csv(path, sep=r"\s+|\t+|,", engine="python", comment="#", header=None)
    df = df.iloc[:, : len(names)]
    df.columns = names
    for c in names:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna().sort_values(names[0]).reset_index(drop=True)


def read_ccd(path: Path) -> pd.DataFrame:
    df = _read2(path, ["Age_Ma", "Global_CCD_m", "CCD_basin_min_m", "CCD_basin_max_m"][:4])
    df["CCD_m_pos_down"] = -df["Global_CCD_m"]
    return df


def read_sea(path: Path) -> pd.DataFrame:
    return _read2(path, ["Age_Ma", "Sea_level_m"])


def match_on_ccd_ages(ccd: pd.DataFrame, sea: pd.DataFrame) -> pd.DataFrame:
    ages = ccd["Age_Ma"].to_numpy(float)
    out = ccd.copy()
    out["Sea_level_m"] = np.interp(ages, sea["Age_Ma"].to_numpy(float), sea["Sea_level_m"].to_numpy(float))
    return out


def pearson_with_ar1_p(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = len(x)
    r, p_naive = stats.pearsonr(x, y)
    if n < 4:
        return r, p_naive, float(n), np.nan
    r1x = stats.pearsonr(x[:-1], x[1:]).statistic
    r1y = stats.pearsonr(y[:-1], y[1:]).statistic
    neff = float(np.clip(n * (1 - r1x * r1y) / (1 + r1x * r1y), 3.0, n))
    t = r * np.sqrt((neff - 2) / max(1e-12, 1 - r * r))
    return float(r), float(p_naive), neff, float(2 * stats.t.sf(abs(t), df=neff - 2))


def fit_ols(x, y):
    s, i, *_ = stats.linregress(x, y)
    return float(s), float(i)


def fit_rma(x, y):
    r = stats.pearsonr(x, y).statistic
    slope = np.sign(r) * np.std(y, ddof=1) / np.std(x, ddof=1)
    return float(slope), float(np.mean(y) - slope * np.mean(x))


def fit_odr(x, y):
    model = odr.Model(lambda beta, xv: beta[0] * xv + beta[1])
    s0, i0 = fit_ols(x, y)
    out = odr.ODR(odr.RealData(x, y), model, beta0=[s0, i0]).run()
    return float(out.beta[0]), float(out.beta[1])


def fit_theilsen(x, y):
    ts = stats.theilslopes(y, x)
    return float(ts.slope), float(ts.intercept)


METHODS: Dict[str, Callable] = {"OLS": fit_ols, "RMA": fit_rma, "ODR": fit_odr, "Theil-Sen_robust": fit_theilsen}


def summarize_fit(method, x, y, fit_func, notes=""):
    slope, intercept = fit_func(x, y)
    r, p_naive, neff, p_ar1 = pearson_with_ar1_p(x, y)
    return FitResult(method, slope, intercept, r, r * r, p_naive, len(x), neff, p_ar1, notes)


def detrend_linear(v, age):
    s, i = fit_ols(age, v)
    return v - (s * age + i)


def moving_block_bootstrap(age, x, y, n_boot=10000, block_len=7, seed=42):
    rng = np.random.default_rng(seed)
    n = len(x)
    starts = np.arange(0, n - block_len + 1)
    n_blocks = int(np.ceil(n / block_len))
    rows = []
    for _ in range(n_boot):
        idx = np.concatenate([np.arange(s, s + block_len) for s in rng.choice(starts, size=n_blocks, replace=True)])[:n]
        xb, yb = x[idx], y[idx]
        if np.std(xb) == 0 or np.std(yb) == 0:
            continue
        rb = stats.pearsonr(xb, yb).statistic
        for name, f in METHODS.items():
            try:
                s, i = f(xb, yb)
                rows.append({"method": name, "slope": s, "intercept": i, "r": rb, "r2": rb * rb})
            except Exception:
                continue
    return pd.DataFrame(rows)


def ci_table(boot):
    rows = []
    for method, g in boot.groupby("method"):
        for var in ["slope", "intercept", "r", "r2"]:
            q = np.nanpercentile(g[var], [2.5, 50, 97.5])
            rows.append({"method": method, "parameter": var, "p2.5": q[0], "median": q[1], "p97.5": q[2]})
    return pd.DataFrame(rows)


def lag_correlation(age, x, y, max_lag_ma=10):
    rows = []
    for lag in range(-max_lag_ma, max_lag_ma + 1):
        xs = np.interp(age, age + lag, x, left=np.nan, right=np.nan)
        m = np.isfinite(xs) & np.isfinite(y)
        if m.sum() >= 8:
            r, _, neff, p_ar1 = pearson_with_ar1_p(xs[m], y[m])
            rows.append({"lag_Ma_sea_level_leads_CCD": lag, "r": r, "r2": r * r, "n": int(m.sum()), "neff_ar1": neff, "p_ar1": p_ar1})
    return pd.DataFrame(rows)


def run(ccd_path, sea_path, outdir, n_boot=10000, block_len=7, max_lag=10) -> Dict:
    ccd_path, sea_path, outdir = Path(ccd_path), Path(sea_path), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = match_on_ccd_ages(read_ccd(ccd_path), read_sea(sea_path))
    df.to_csv(outdir / "matched_timeseries.csv", index=False)
    age = df["Age_Ma"].to_numpy(float)
    x = df["Sea_level_m"].to_numpy(float)
    y = df["CCD_m_pos_down"].to_numpy(float)

    raw = pd.DataFrame([summarize_fit(n, x, y, f, "raw").__dict__ for n, f in METHODS.items()])
    raw.to_csv(outdir / "regression_results.csv", index=False)

    x_dt, y_dt = detrend_linear(x, age), detrend_linear(y, age)
    det = pd.DataFrame([summarize_fit(n, x_dt, y_dt, f, "detrended").__dict__ for n, f in METHODS.items()])
    det.to_csv(outdir / "detrended_regression_results.csv", index=False)

    dx, dy = np.diff(x), np.diff(y)
    dif = pd.DataFrame([summarize_fit(n, dx, dy, f, "first-difference").__dict__ for n, f in METHODS.items()])
    dif.to_csv(outdir / "first_difference_regression_results.csv", index=False)

    boot = moving_block_bootstrap(age, x, y, n_boot=n_boot, block_len=block_len)
    boot_ci = ci_table(boot)
    boot_ci.to_csv(outdir / "bootstrap_95CI.csv", index=False)

    # Pointwise 95% band for the preferred RMA calibration, built from the JOINT
    # replicates. The marginal percentiles in bootstrap_95CI.csv cannot give this: slope
    # and intercept are strongly anti-correlated in a fit, so pairing their extremes
    # would draw a band several times too wide. Each replicate is evaluated as a line
    # over the observed sea-level range and the percentiles taken down the columns.
    _rma_boot = boot[boot["method"] == "RMA"]
    _xg = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    _lines = (_rma_boot["slope"].to_numpy()[:, None] * _xg[None, :]
              + _rma_boot["intercept"].to_numpy()[:, None])
    _q = np.nanpercentile(_lines, [2.5, 50, 97.5], axis=0)
    pd.DataFrame({"Sea_level_m": _xg, "CCD_p2.5": _q[0],
                  "CCD_median": _q[1], "CCD_p97.5": _q[2]}).to_csv(
        outdir / "bootstrap_band_rma.csv", index=False)

    lag = lag_correlation(age, x, y, max_lag_ma=max_lag)
    lag.to_csv(outdir / "lag_correlation.csv", index=False)

    with open(outdir / "regression_summary.txt", "w") as f:
        f.write("Sea level vs CCD regression summary (CCD positive-down)\n\n")
        f.write("RAW\n" + raw.to_string(index=False) + "\n\n")
        f.write("DETRENDED\n" + det.to_string(index=False) + "\n\n")
        f.write("FIRST DIFFERENCE\n" + dif.to_string(index=False) + "\n\n")
        f.write("BOOTSTRAP 95% CI (raw)\n" + boot_ci.to_string(index=False) + "\n")

    from ccdworkflow.io import save_matplotlib_figure

    # preferred RMA scatter figure
    rma = raw.loc[raw["method"] == "RMA"].iloc[0]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(x, y, s=40, facecolors="white", edgecolors="black", zorder=2, label="Matched 1 Myr values")
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 200)
    ax.plot(xx, rma["slope"] * xx + rma["intercept"], lw=2, color="black", label="RMA fit")
    ax.invert_yaxis()
    ax.set_xlabel("Sea level (m)"); ax.set_ylabel("CCD (m; deeper downward)")
    ax.set_title("Preferred RMA regression, 0-52 Ma")
    ax.text(0.03, 0.97, f"slope={rma['slope']:.2f} m/m\nintercept={rma['intercept']:.0f} m\n"
            f"r={rma['r']:.3f}, R2={rma['r2']:.3f}", transform=ax.transAxes, va="top",
            bbox=dict(boxstyle="round", fc="white", ec="black", alpha=0.9))
    ax.legend(frameon=False)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step6_preferred_rma_regression", dpi=300, step="step6")
    plt.close(fig)

    # sensitivity of the fit to regression method (all four fits over the data)
    _styles = {"OLS": ("#1f77b4", "-"), "RMA": ("#000000", "-"),
               "ODR": ("#d62728", "--"), "Theil-Sen_robust": ("#2ca02c", ":")}
    _labels = {"OLS": "OLS", "RMA": "RMA (preferred)", "ODR": "ODR", "Theil-Sen_robust": "Theil-Sen"}
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    ax.scatter(x, y, s=40, facecolors="none", edgecolors="0.25", linewidths=1.1, zorder=3, label="Matched 1 Myr values")
    for _, rw in raw.iterrows():
        c, ls = _styles[rw["method"]]
        lw = 2.8 if rw["method"] == "RMA" else 1.7
        ax.plot(xx, rw["slope"] * xx + rw["intercept"], color=c, ls=ls, lw=lw,
                label=f"{_labels[rw['method']]}: {rw['slope']:.2f} m/m, {rw['intercept']:.0f} m")
    ax.invert_yaxis()
    ax.set_xlabel("Sea level (m)"); ax.set_ylabel("CCD (m; deeper downward)")
    ax.set_title(f"Sensitivity to regression method, 0-52 Ma (r={rma['r']:.2f}, R2={rma['r2']:.2f})")
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    # This diagnostic IS Figure S1 of the supplement, so it goes straight into the
    # paper's figure folder as well - no manual copy step.
    save_matplotlib_figure(fig, "step6_regression_method_sensitivity", dpi=300, step="step6",
                           paper_name="FigS1_regression_method_sensitivity")
    plt.close(fig)

    # detrended scatter (shared long-term trend removed from both records)
    rma_dt = det.loc[det["method"] == "RMA"].iloc[0]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(x_dt, y_dt, s=36, facecolors="white", edgecolors="black", zorder=2, label="Detrended values")
    xxd = np.linspace(np.nanmin(x_dt), np.nanmax(x_dt), 200)
    ax.plot(xxd, rma_dt["slope"] * xxd + rma_dt["intercept"], lw=2, color="black", label="RMA fit")
    ax.set_xlabel("Detrended sea level (m)"); ax.set_ylabel("Detrended CCD (m)")
    ax.set_title(f"Detrended sea level vs CCD, 0-52 Ma (r={rma_dt['r']:.2f}, R2={rma_dt['r2']:.2f}, p_AR1={rma_dt['p_ar1']:.2f})")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step6_detrended_scatter", dpi=300, step="step6")
    plt.close(fig)

    # first-difference scatter (1 Myr changes)
    rma_df = dif.loc[dif["method"] == "RMA"].iloc[0]
    fig, ax = plt.subplots(figsize=(6.5, 5))
    ax.scatter(dx, dy, s=36, facecolors="white", edgecolors="black", zorder=2, label="1 Myr changes")
    xxf = np.linspace(np.nanmin(dx), np.nanmax(dx), 200)
    ax.plot(xxf, rma_df["slope"] * xxf + rma_df["intercept"], lw=2, color="black", label="RMA fit")
    ax.set_xlabel("Δ sea level (m / Myr)"); ax.set_ylabel("Δ CCD (m / Myr)")
    ax.set_title(f"First-difference sea level vs CCD, 0-52 Ma (r={rma_df['r']:.2f}, R2={rma_df['r2']:.2f})")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step6_first_difference_scatter", dpi=300, step="step6")
    plt.close(fig)

    # lag-correlation curve
    if len(lag):
        fig, ax = plt.subplots(figsize=(7.0, 4.2))
        ax.axhline(0, color="0.7", lw=0.8); ax.axvline(0, color="0.7", lw=0.8)
        ax.plot(lag["lag_Ma_sea_level_leads_CCD"], lag["r"], "-o", color="#b3202c", ms=4)
        ax.set_xlabel("Lag (Myr; + = sea level leads CCD)"); ax.set_ylabel("Correlation r")
        ax.set_title("Lag correlation between long-term sea level and CCD, 0-52 Ma")
        fig.tight_layout()
        save_matplotlib_figure(fig, "step6_lag_correlation", dpi=300, step="step6")
        plt.close(fig)

    return {"raw": raw, "rma_slope": float(rma["slope"]), "rma_intercept": float(rma["intercept"])}
