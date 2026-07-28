#!/usr/bin/env python3
"""
Step 5 — Sea-level ↔ CCD regression and pre-Cenozoic CCD prediction.

Two parts:

1. **Regression diagnostics** (0–52 Ma) via :mod:`regression_diagnostics` — OLS /
   RMA (preferred) / ODR / Theil–Sen, detrended and first-difference tests,
   AR(1) p-values, moving-block bootstrap, lag test. Written to
   ``outputs/step6_regression/diagnostics/``.

2. **Prediction** of CCD over 0–205 Ma from the sea-level envelope, using the
   FIXED preferred calibration ``Global_CCD = slope*SL + intercept``
   (slope = 5.72 m/m, intercept = -4456.8 m). A 95% *trend* uncertainty envelope
   is derived from the residual RMSE over the 0–52 Ma calibration interval,
   scaled by the effective degrees of freedom (Neff = cal_window / Lcorr).
   Reproduces ``5_regression/2_predict_CCD_regress.sh`` but in pure Python
   (Akima resampling instead of ``gmt sample1d``).

Inputs:
    outputs/step5_ccd_lowpass/low_pass_filtered_ccd_0-52Ma.txt
    outputs/step4_sealevel_envelope/sea_level_quantile_envelope_0-205Ma.txt  (prediction)
    outputs/step4_sealevel_envelope/sea_level_quantile_envelope_0-52Ma.txt   (regression)

Output:
    outputs/step6_regression/predicted_ccd_sl_0-205Ma.txt
        Age_Ma, Predicted_CCD_m, Lo95_m, Hi95_m
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import Akima1DInterpolator

from ccdworkflow import config
from ccdworkflow.io import read_table, read_xy, write_table

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # sibling module in this step folder
import regression_diagnostics as rdiag

P = config.PARAMS


def _akima_grid(age, val, tstart, tend, dt):
    age = np.asarray(age, float)
    val = np.asarray(val, float)
    order = np.argsort(age)
    age, val = age[order], val[order]
    keep = np.concatenate([[True], np.diff(age) > 0])   # strictly increasing
    age, val = age[keep], val[keep]
    grid = np.arange(tstart, tend + dt, dt)
    interp = Akima1DInterpolator(age, val)
    return grid, interp(grid)


def predict(slope: float, intercept: float) -> pd.DataFrame:
    sl = read_xy(config.SL_ENVELOPE_FULL, names=["Age_Ma", "SL_m"])
    grid, sl_grid = _akima_grid(sl["Age_Ma"], sl["SL_m"], P.pred_tstart, P.pred_tend, P.pred_dt)

    # predicted central CCD (negative-up Global_CCD convention)
    pred = slope * sl_grid + intercept

    # residual RMSE over calibration interval vs observed low-pass CCD
    obs = read_table(config.LOWPASS_CCD_0_52, names=["Age_Ma", "Global_CCD_m", "lo", "hi"])
    obs = obs.dropna(subset=["Age_Ma", "Global_CCD_m"])
    _, obs_grid = _akima_grid(obs["Age_Ma"], obs["Global_CCD_m"], P.pred_tstart, P.pred_tend, P.pred_dt)
    cal = (grid >= 0.0) & (grid <= P.ccd_cal_max_age) & np.isfinite(obs_grid)
    resid = pred[cal] - obs_grid[cal]
    rmse = float(np.sqrt(np.mean(resid ** 2)))
    neff = (P.ccd_cal_max_age - 0.0) / P.pred_corr_lengthscale_myr
    sigma_trend = rmse / np.sqrt(neff)
    half = P.pred_z95 * sigma_trend
    print(f"[step6] prediction RMSE={rmse:.2f} m  Neff={neff:.2f}  sigma_trend={sigma_trend:.2f} m  95%=±{half:.1f} m")

    return pd.DataFrame({
        "Age_Ma": grid,
        "Predicted_CCD_m": pred,
        "Lo95_m": pred - half,
        "Hi95_m": pred + half,
    })


def main() -> None:
    config.ensure_dirs()

    # 1. regression diagnostics
    diag_dir = config.STEP6_DIR / "diagnostics"
    res = rdiag.run(config.LOWPASS_CCD_0_52, config.SL_ENVELOPE_0_52, diag_dir)
    print(f"[step6] regression diagnostics -> {diag_dir}  "
          f"(RMA slope={res['rma_slope']:.3f}, intercept={res['rma_intercept']:.1f}; "
          f"CCD positive-down)")

    # 2. prediction 0-205 Ma. Convert positive-down RMA fit to the negative-up
    #    Global_CCD convention used for prediction: Global_CCD = -CCD.
    if P.use_fitted_rma:
        slope = -res["rma_slope"]
        intercept = -res["rma_intercept"]
        print(f"[step6] using fitted RMA -> pred_slope={slope:.3f}, pred_intercept={intercept:.1f}")
    else:
        slope, intercept = P.pred_slope, P.pred_intercept
        print(f"[step6] using published fixed pred_slope={slope}, pred_intercept={intercept}")
    pred_df = predict(slope, intercept)
    out = write_table(config.PREDICTED_CCD_SL, pred_df,
                      header=["Age_Ma", "Predicted_CCD_m", "Lo95_m", "Hi95_m"],
                      float_format="%.3f")
    print(f"[step6] wrote {out}  ({len(pred_df)} rows)")


if __name__ == "__main__":
    main()
