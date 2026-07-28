#!/usr/bin/env python3
"""
Step 0 — Regional CCD emendation / resampling.

Resample the raw 0.5 Myr regional CCD curves (Atlantic, Pacific, Indian) onto a
uniform 1 Myr grid using Akima interpolation, matching the original
``resample_{ATL,IND,PAC}_CCD_1my.py`` scripts (Akima, integer output).

Inputs  (data/regional_ccd/):
    ATL_CCD_0.5my.tsv, PAC_CCD_0.5my.tsv, IND_CCD_0.5my.tsv
        two columns: Age_Ma, Mean_CCD_m  (commented header)

Outputs (outputs/step1_regional_ccd/):
    ATL_CCD_1my_akima.txt, PAC_CCD_1my_akima.txt, IND_CCD_1my_akima.txt

Note
----
These resampled regional curves, together with the ocean-basin area fractions
from step 1, are assembled (by hand, in ``ccds_and_oceanbasin_fractions.xlsx``)
into the input read by step 2. This step is included for completeness and
reproducibility of the regional curves; the xlsx remains the manual join point.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import Akima1DInterpolator

from ccdworkflow import config
from ccdworkflow.io import read_xy, write_table

BASINS = {
    "ATL": config.REGIONAL_CCD / "ATL_CCD_0.5my.tsv",
    "PAC": config.REGIONAL_CCD / "PAC_CCD_0.5my.tsv",
    "IND": config.REGIONAL_CCD / "IND_CCD_0.5my.tsv",
}


def resample_basin(src, dt: float = 1.0):
    df = read_xy(src, names=["Age_Ma", "CCD_m"])
    age = df["Age_Ma"].to_numpy(float)
    ccd = df["CCD_m"].to_numpy(float)
    new_age = np.arange(np.ceil(age.min()), np.floor(age.max()) + dt, dt)
    interp = Akima1DInterpolator(age, ccd)
    new_ccd = np.rint(interp(new_age)).astype(int)
    return np.column_stack([new_age.astype(int), new_ccd])


def _plot(collected: dict) -> None:
    """Diagnostic figure: raw 0.5 Myr regional CCDs vs the 1 Myr Akima resample."""
    if not collected:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ccdworkflow.io import save_matplotlib_figure

    colours = {"ATL": "#1f77b4", "PAC": "#d62728", "IND": "#2ca02c"}
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for name, (raw, res) in collected.items():
        c = colours.get(name, None)
        ax.plot(raw["Age_Ma"], raw["CCD_m"], color=c, lw=0.7, alpha=0.4)
        ax.plot(res[:, 0], res[:, 1], color=c, lw=1.8, label=f"{name} (1 Myr Akima)")
    ax.invert_xaxis()
    ax.set_xlabel("Age (Ma)"); ax.set_ylabel("Regional CCD (m)")
    ax.set_title("Regional CCD curves: raw 0.5 Myr (faint) and 1 Myr Akima resample")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_matplotlib_figure(fig, "step1_regional_ccd_resample", dpi=300, step="step1")
    plt.close(fig)


def main() -> None:
    config.ensure_dirs()
    import pandas as pd
    from ccdworkflow.io import read_xy

    collected: dict = {}
    for name, src in BASINS.items():
        if not src.exists():
            print(f"[step1] SKIP {name}: missing input {src}")
            continue
        arr = resample_basin(src)
        df = pd.DataFrame({"Age_Ma": arr[:, 0], "CCD_m": arr[:, 1]})
        out = write_table(config.STEP1_DIR / f"{name}_CCD_1my_akima.txt", df,
                          header=["Age_Ma", "CCD_m"], float_format="%.0f")
        print(f"[step1] wrote {out}")
        collected[name] = (read_xy(src, names=["Age_Ma", "CCD_m"]), arr)
    _plot(collected)


if __name__ == "__main__":
    main()
