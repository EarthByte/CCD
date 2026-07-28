#!/usr/bin/env python3
"""
Step 6 — Hybrid CCD: splice observed + sea-level-predicted CCD, add planktogenic term.

Build the final hybrid CCD curve:

1. Read observed low-pass CCD (0–52 Ma) and sea-level-predicted CCD (0–205 Ma).
2. Add a planktogenic correction to the predicted curve (mean and both bounds):
   ``corr = slope * (Age - start_age)`` for ``Age >= start_age``
   (slope = 20 m/Myr, start_age = 115 Ma — see README "Known discrepancies").
3. Splice: observed for Age <= 52, a linear bridge node at 53, predicted for Age >= 54.
4. Write the combined curve and a closed uncertainty-envelope polygon.

Inputs:
    outputs/step5_ccd_lowpass/low_pass_filtered_ccd_0-52Ma.txt
    outputs/step6_regression/predicted_ccd_sl_0-205Ma.txt

Outputs (outputs/step7_hybrid/):
    hybrid_ccd_obs_pred_combined.txt   Age_Ma, CCD_m, CCD_min_m, CCD_max_m
    hybrid_ccd_error_envelope.txt      closed polygon for shaded plotting
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ccdworkflow import config
from ccdworkflow.io import read_table, write_table, write_polygon

P = config.PARAMS


def _order_bounds(df, lo_col, hi_col):
    lo = np.minimum(df[lo_col], df[hi_col])
    hi = np.maximum(df[lo_col], df[hi_col])
    return lo, hi


def read_observed():
    df = read_table(config.LOWPASS_CCD_0_52,
                    names=["Age", "Mean", "b1", "b2"]).dropna().sort_values("Age")
    lo, hi = _order_bounds(df, "b1", "b2")
    return pd.DataFrame({"Age": df["Age"].to_numpy(), "Mean": df["Mean"].to_numpy(),
                         "Min": lo.to_numpy(), "Max": hi.to_numpy()})


def read_predicted():
    df = read_table(config.PREDICTED_CCD_SL,
                    names=["Age", "Mean", "b1", "b2"]).dropna().sort_values("Age")
    lo, hi = _order_bounds(df, "b1", "b2")
    return pd.DataFrame({"Age": df["Age"].to_numpy(), "Mean": df["Mean"].to_numpy(),
                         "Min": lo.to_numpy(), "Max": hi.to_numpy()})


def apply_planktogenic(df):
    out = df.copy()
    age = out["Age"].to_numpy()
    corr = np.where(age >= P.plankton_start_age_ma,
                    P.plankton_slope_m_per_myr * (age - P.plankton_start_age_ma), 0.0)
    for c in ("Mean", "Min", "Max"):
        out[c] = out[c] + corr
    return out


def build_hybrid(obs, pred):
    obs_part = obs[obs["Age"] <= P.splice_obs_max_age]
    obs52 = obs.loc[np.isclose(obs["Age"], P.splice_obs_max_age)]
    pred54 = pred.loc[np.isclose(pred["Age"], P.splice_pred_min_age)]
    if obs52.empty or pred54.empty:
        raise ValueError("Need observed at 52 Ma and predicted at 54 Ma for the bridge.")
    o, q = obs52.iloc[0], pred54.iloc[0]
    f = 0.5
    bridge = pd.DataFrame([{
        "Age": P.splice_bridge_age,
        "Mean": o["Mean"] + f * (q["Mean"] - o["Mean"]),
        "Min": o["Min"] + f * (q["Min"] - o["Min"]),
        "Max": o["Max"] + f * (q["Max"] - o["Max"]),
    }])
    pred_part = pred[pred["Age"] >= P.splice_pred_min_age]
    return pd.concat([obs_part, bridge, pred_part], ignore_index=True).sort_values("Age").reset_index(drop=True)


def main() -> None:
    config.ensure_dirs()
    obs = read_observed()
    pred = apply_planktogenic(read_predicted())
    hyb = build_hybrid(obs, pred)

    out = pd.DataFrame({
        "Age_Ma": hyb["Age"],
        "CCD_m": hyb["Mean"],
        "CCD_min_m": hyb["Min"],
        "CCD_max_m": hyb["Max"],
    })
    path = write_table(config.HYBRID_CCD, out,
                       header=["Age_Ma", "CCD_m", "CCD_min_m", "CCD_max_m"],
                       float_format="%.2f")
    env = write_polygon(config.HYBRID_CCD_ENVELOPE,
                        hyb["Age"].to_numpy(), hyb["Min"].to_numpy(), hyb["Max"].to_numpy())
    print(f"[step7] wrote {path}  ({len(out)} rows)")
    print(f"[step7] wrote {env}")

    # diagnostic figure: hybrid CCD with uncertainty + its observed/predicted parts
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ccdworkflow.io import save_matplotlib_figure
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.fill_between(hyb["Age"], hyb["Min"], hyb["Max"], color="0.75", alpha=0.6, label="uncertainty")
    ax.plot(hyb["Age"], hyb["Mean"], color="#b3202c", lw=2.4, label="Hybrid CCD")
    ax.plot(obs["Age"], obs["Mean"], color="#2ca02c", lw=1.6, label="Observed (≤52 Ma)")
    pred_part = pred[pred["Age"] >= P.splice_pred_min_age]
    ax.plot(pred_part["Age"], pred_part["Mean"], color="#3f74ad", lw=1.6, label="Sea-level prediction (≥54 Ma)")
    ax.invert_xaxis()
    ax.set_xlabel("Age (Ma)"); ax.set_ylabel("CCD (m)")
    ax.set_title("Hybrid CCD reconstruction with uncertainty, 0–205 Ma")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step7_hybrid_ccd_components", dpi=300, step="step7")
    plt.close(fig)


if __name__ == "__main__":
    main()
