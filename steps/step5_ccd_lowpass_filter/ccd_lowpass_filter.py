#!/usr/bin/env python3
"""
Step 4 — Low-pass filter the global CCD to match the sea-level bandwidth.

Apply a Butterworth low-pass filter (zero-phase, ``filtfilt``) to the global CCD
central curve and to its inter-basin dispersion bounds, so that the CCD and the
long-term sea-level envelope are compared at a similar effective bandwidth.

Filter parameters (config.PARAMS), tuned in the original interactive notebook and
fixed here: cutoff = 0.5, order = 1, no post-smoothing (moving-average width 1),
fs = 1 sample/Myr. The cutoff is normalised by ``nyquist = fs`` (i.e. Wn = 0.5),
reproducing the original ``butter_lowpass`` exactly.

Input:
    outputs/step3_global_ccd/global_ccd_with_basin_dispersion_0-52Ma.txt

Output:
    outputs/step5_ccd_lowpass/low_pass_filtered_ccd_0-52Ma.txt
        Age_Ma, Global_CCD_m, CCD_minus_dispersion_m, CCD_plus_dispersion_m
        (age integer; values to 1 decimal; commented header)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

from ccdworkflow import config
from ccdworkflow.io import read_table, write_table

P = config.PARAMS


def _butter_lowpass(cutoff: float, fs: float, order: int):
    nyquist = 1.0 * fs           # NB: original uses fs (not fs/2) as the divisor
    wn = min(cutoff / nyquist, 0.99)
    return butter(order, wn, btype="low", analog=False)


def lowpass(signal: np.ndarray) -> np.ndarray:
    b, a = _butter_lowpass(P.ccd_butter_cutoff, P.ccd_fs, P.ccd_butter_order)
    out = filtfilt(b, a, signal)
    if P.ccd_smooth_samples > 1:
        k = P.ccd_smooth_samples
        out = np.convolve(out, np.ones(k) / k, mode="same")
    return out


def main() -> None:
    config.ensure_dirs()
    df = read_table(config.GLOBAL_CCD_0_52,
                    names=["Age_Ma", "Global_CCD_m",
                           "CCD_minus_dispersion_m", "CCD_plus_dispersion_m"])
    df = df.dropna(subset=["Age_Ma", "Global_CCD_m"]).reset_index(drop=True)

    age = df["Age_Ma"].to_numpy(float)
    central = lowpass(df["Global_CCD_m"].to_numpy(float))
    minus = lowpass(df["CCD_minus_dispersion_m"].to_numpy(float))
    plus = lowpass(df["CCD_plus_dispersion_m"].to_numpy(float))

    out = pd.DataFrame({
        "Age_Ma": np.rint(age).astype(int),
        "Global_CCD_m": np.round(central, 1),
        "CCD_minus_dispersion_m": np.round(minus, 1),
        "CCD_plus_dispersion_m": np.round(plus, 1),
    })
    path = write_table(config.LOWPASS_CCD_0_52, out,
                       header=["Age_Ma", "Global_CCD_m",
                               "CCD_minus_dispersion_m", "CCD_plus_dispersion_m"],
                       float_format="%.1f")
    print(f"[step5] wrote {path}  ({len(out)} rows)")

    # diagnostic figure: raw vs low-pass-filtered global CCD
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ccdworkflow.io import save_matplotlib_figure
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.fill_between(age, minus, plus, color="0.6", alpha=0.3, label="filtered dispersion")
    ax.plot(age, df["Global_CCD_m"].to_numpy(float), color="#9aa7b4", lw=0.9, alpha=0.9, label="raw global CCD")
    ax.plot(age, central, color="black", lw=2.4, label="low-pass filtered CCD")
    ax.invert_xaxis()
    ax.set_xlabel("Age (Ma)"); ax.set_ylabel("Global CCD (m)")
    ax.set_title(f"Butterworth low-pass CCD (cutoff={P.ccd_butter_cutoff}, order={P.ccd_butter_order}), 0–52 Ma")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step5_raw_vs_lowpass_ccd", dpi=300, step="step5")
    plt.close(fig)


if __name__ == "__main__":
    main()
