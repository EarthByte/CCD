#!/usr/bin/env python3
"""
Standalone reproduction of atmosphere_influx_comparison_sameaxes.png from
Notebook 05, one figure per CCD model.

Eight series, all from 05_atmospheric_influx_all_sources.csv:
    Mid-ocean ridges       — ridge_outflux_*
    Subduction zones       — subduction_outflux_*
    Carbonate platforms    — carbonate_platform_outflux_*
    Rifts (corrected)      — rift_outflux_non_biased_*
    Intraplate volcanism   — intraplate_volcanism_outflux_*
    Total (without rifts)  — gross_atmospheric_outflux_no_rift_*
    Total (with rifts)     — gross_atmospheric_outflux_unbiased_rift_*

Source: cell 30 of 05-Atmospheric-Carbon.ipynb.

Outputs (under ./output/):
    atmosphere_influx_comparison_DM2026.png/pdf
    atmosphere_influx_comparison_BW1991.png/pdf
"""
from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import _common as cfg


COMPONENTS = [
    # CSV stem,                       label,                 colour
    ("ridge_outflux",                 "Mid-ocean ridges",    "dodgerblue"),
    ("subduction_outflux",            "Subduction zones",    "C1"),
    ("carbonate_platform_outflux",    "Carbonate platforms", "C2"),
    ("rift_outflux_non_biased",       "Rifts",               "C3"),
    ("intraplate_volcanism_outflux",  "Intraplate volcanism","C8"),
]


def plot_one(model: str):
    df = cfg.load_atmospheric_influx(model)
    times = df.index.to_numpy()
    max_time = float(times.max())

    fig, ax = plt.subplots(figsize=(8.0, 4.5))

    for stem, label, colour in COMPONENTS:
        ax.fill_between(times,
                        df[f"{stem}_min"].to_numpy(),
                        df[f"{stem}_max"].to_numpy(),
                        color=colour, alpha=0.5)
        ax.plot(times, df[f"{stem}_mean"].to_numpy(),
                color=colour, label=label)

    # Total (without rifts) — black line over grey shading.
    ax.fill_between(times,
                    df["gross_atmospheric_outflux_no_rift_min"].to_numpy(),
                    df["gross_atmospheric_outflux_no_rift_max"].to_numpy(),
                    color="0.65", alpha=0.5)
    ax.plot(times, df["gross_atmospheric_outflux_no_rift_mean"].to_numpy(),
            color="k", label="Total (without rifts)")

    # Total (with rifts) — purple over light purple.
    ax.fill_between(times,
                    df["gross_atmospheric_outflux_unbiased_rift_min"].to_numpy(),
                    df["gross_atmospheric_outflux_unbiased_rift_max"].to_numpy(),
                    color="C4", alpha=0.5)
    ax.plot(times, df["gross_atmospheric_outflux_unbiased_rift_mean"].to_numpy(),
            color="C4", label="Total (with rifts)")

    ax.set_xlim(max_time, 0)
    ax.set_xlabel("Age (Ma)")
    ax.set_ylabel("Atmosphere influx (Mt C/yr)")
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.set_ticks(np.arange(0., max_time + 1, 20.), minor=False)
    ax.tick_params(axis="both", direction="in", which="both",
                   right=True, top=True)
    ax.grid(alpha=0.1, which="both")
    ax.legend(frameon=False, loc="upper left",
              prop={"size": 8}, ncols=2)
    ax.set_title(f"Atmosphere influx — {model}", fontsize=11)

    cfg.save_figure(fig, f"atmosphere_influx_comparison_{model}")
    plt.close(fig)


def main():
    for model in cfg.MODELS:
        print(f"[{model}]")
        plot_one(model)


if __name__ == "__main__":
    main()
