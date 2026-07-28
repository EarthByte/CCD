#!/usr/bin/env python3
"""
Step 3 — Long-term sea-level curve from a rolling high-quantile peak envelope.

Apply a rolling high-quantile, peak-following filter to the short-term sea-level
compilation to obtain a long-term, low-pass, peak-following envelope. This is the
canonical "quantile envelope" method; the experimental peak-spline (PchipInterpolator
+ Savitzky–Golay) alternative from the original notebook has been REMOVED, as only
the quantile envelope was used.

Parameters (config.PARAMS): window = 2.0 Myr, quantile = 0.90, smoothing = 1.5 Myr.

Input:
    data/sealevel/Miller_Haq_SeaLevel_ShortTerm_hybrid.tsv   (Age_Ma, SL_m; ~0.1 Myr)

Outputs (outputs/step4_sealevel_envelope/):
    sea_level_quantile_envelope_0-205Ma.txt   (full range; feeds step 5 prediction)
    sea_level_quantile_envelope_0-52Ma.txt    (calibration subset; feeds steps 4 & 5)

Optional validation (``--validate``) reproduces the comparison against the
published Haq long-term curve, reporting Pearson r, RMSE, MAE and bias.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from ccdworkflow import config
from ccdworkflow.io import read_xy, write_table

P = config.PARAMS
RAW_SL = config.SEALEVEL / "Miller_Haq_SeaLevel_ShortTerm_hybrid.tsv"   # envelope input
HAQ_SHORTTERM = config.SEALEVEL / "Haq_SeaLevel_ShortTerm_hybrid.tsv"   # validation input
HAQ_LONGTERM = config.SEALEVEL / "Haq87_Longterm_v3.txt"                # validation target


def _to_samples(myr: float, dt: float) -> int:
    n = int(round(max(myr, dt) / dt))
    if n % 2 == 0:
        n += 1
    return max(n, 3)


def quantile_envelope(ages: np.ndarray, sl: np.ndarray,
                      window_myr: float, quantile: float, smooth_myr: float) -> np.ndarray:
    """Rolling high-quantile peak-following envelope (canonical method)."""
    diffs = np.diff(ages)
    dt = np.median(diffs[diffs != 0]) if np.any(diffs != 0) else 1.0
    win_q = _to_samples(window_myr, dt)
    s = pd.Series(sl, dtype=float)
    env = s.rolling(win_q, center=True, min_periods=max(3, win_q // 4)).quantile(quantile)
    env = env.interpolate(limit_direction="both").bfill().ffill()
    if smooth_myr > 0:
        win_s = _to_samples(smooth_myr, dt)
        env = env.rolling(win_s, center=True, min_periods=1).mean().bfill().ffill()
    return env.to_numpy()


def compute() -> pd.DataFrame:
    df = read_xy(RAW_SL, names=["Age_Ma", "SL_m"])
    ages = df["Age_Ma"].to_numpy(float)
    sl = df["SL_m"].to_numpy(float)
    env = quantile_envelope(ages, sl, P.sl_window_myr, P.sl_quantile, P.sl_smooth_myr)
    return pd.DataFrame({"Age_Ma": ages, "SL_envelope_m": env})


def _plot_envelope(env: pd.DataFrame) -> None:
    """Diagnostic figure: short-term sea-level record with the peak-following
    quantile envelope used downstream."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ccdworkflow.io import save_matplotlib_figure

    raw = read_xy(RAW_SL, names=["Age_Ma", "SL_m"])
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    ax.plot(raw["Age_Ma"], raw["SL_m"], color="#9aa7b4", lw=0.7, alpha=0.8, label="Short-term sea level (input)")
    ax.plot(env["Age_Ma"], env["SL_envelope_m"], color="#b3202c", lw=2.2,
            label=f"0.{int(P.sl_quantile*100)} quantile peak-following envelope")
    ax.set_xlim(env["Age_Ma"].max(), env["Age_Ma"].min())
    ax.set_xlabel("Age (Ma)"); ax.set_ylabel("Sea level (m)")
    ax.set_title("Long-term sea-level envelope (rolling high-quantile, peak-following)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step4_sealevel_quantile_envelope", dpi=300, step="step4")
    plt.close(fig)


def validate() -> None:
    # Validation applies the operator to the Haq short-term record and compares
    # against the published Haq long-term curve (reproduces the notebook check).
    df = read_xy(HAQ_SHORTTERM, names=["Age_Ma", "SL_m"])
    haq = read_xy(HAQ_LONGTERM, names=["Age_Ma", "SL_m"])
    a0, a1 = df["Age_Ma"].min(), df["Age_Ma"].max()
    grid = np.arange(a0, a1 + P.sl_grid_dt, P.sl_grid_dt)
    sl_u = np.interp(grid, df["Age_Ma"], df["SL_m"])
    env = quantile_envelope(grid, sl_u, P.sl_window_myr, P.sl_quantile, P.sl_smooth_myr)
    long_u = np.interp(grid, haq["Age_Ma"], haq["SL_m"])
    lo = max(a0, haq["Age_Ma"].min()) + 2.0
    hi = min(a1, haq["Age_Ma"].max()) - 2.0
    m = (grid >= lo) & (grid <= hi)
    e, l = env[m], long_u[m]
    r = float(np.corrcoef(e, l)[0, 1])
    rmse = float(np.sqrt(np.mean((e - l) ** 2)))
    mae = float(np.mean(np.abs(e - l)))
    bias = float(np.mean(e - l))
    print(f"[step4][validate] vs Haq long-term (edges trimmed): "
          f"r={r:.3f}  RMSE={rmse:.1f} m  MAE={mae:.1f} m  bias={bias:+.1f} m")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ccdworkflow.io import save_matplotlib_figure
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    ax.plot(grid, long_u, color="#3f74ad", lw=2.0, label="Published Haq long-term")
    ax.plot(grid, env, color="#b3202c", lw=1.8, label="Quantile envelope of Haq short-term")
    ax.axvspan(grid.min(), lo, color="0.9"); ax.axvspan(hi, grid.max(), color="0.9")
    ax.set_xlim(grid.max(), grid.min())
    ax.set_xlabel("Age (Ma)"); ax.set_ylabel("Sea level (m)")
    ax.set_title(f"Validation of the peak-following operator (r={r:.2f}, RMSE={rmse:.1f} m, MAE={mae:.1f} m)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step4_validation_vs_haq_longterm", dpi=300, step="step4")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true", help="report fit vs Haq long-term")
    args = ap.parse_args()

    config.ensure_dirs()
    env = compute()
    full = write_table(config.SL_ENVELOPE_FULL, env,
                       header=["Age_Ma", "SL_envelope_m"], float_format="%.3f")
    sub = env[env["Age_Ma"] <= P.ccd_cal_max_age].reset_index(drop=True)
    sub_out = write_table(config.SL_ENVELOPE_0_52, sub,
                          header=["Age_Ma", "SL_envelope_m"], float_format="%.3f")
    print(f"[step4] wrote {full}  ({len(env)} rows)")
    print(f"[step4] wrote {sub_out}  ({len(sub)} rows, 0-52 Ma)")
    _plot_envelope(env)
    if args.validate:
        validate()


if __name__ == "__main__":
    main()
