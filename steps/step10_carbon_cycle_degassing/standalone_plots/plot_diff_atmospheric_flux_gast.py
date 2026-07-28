#!/usr/bin/env python3
"""
Standalone reproduction of 170-0_diff_atmospheric_flux_both_rift_options_gast.png
from Notebook 05, one figure per CCD model.

Two stacked panels (same x-axis = age 170 Ma -> 0 Ma):
    Upper:  Net carbon outflux (corrected / uncorrected rift length) +
            Total outflux       (corrected / uncorrected rift length)
            All four series read from 05_carbon_flux.csv:
                net_carbon_outflux_no_sed_unbiased_rifts   (corrected)
                net_carbon_outflux_no_sed_biased_rifts     (uncorrected)
                gross_atmospheric_influx_no_sed_unbiased_rifts (corrected total)
                gross_atmospheric_influx_no_sed_biased_rifts   (uncorrected total)

    Lower:  Global temperature — Mills GAST from 05_global_temperature.csv
            (grey shaded uncertainty + black mean) over which the Scotese
            curve from utils/GAST_800Ma_Scotese.xlsx is overlaid in blue.

Source: cell 75 of 05-Atmospheric-Carbon.ipynb (j == 0 branch).

Outputs (under ./output/):
    170-0_diff_atmospheric_flux_gast_DM2026.png/pdf
    170-0_diff_atmospheric_flux_gast_BW1991.png/pdf
"""
from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt

import _common as cfg


def _make_panels():
    fig, (ax_flux, ax_gast) = plt.subplots(
        ncols=1, nrows=2, figsize=(10, 7),
        constrained_layout=True, height_ratios=[1, 1],
    )
    return fig, ax_flux, ax_gast


def _fill_and_plot(ax, x, mid, lo, hi, *, color, label, linestyle="-",
                   fill_alpha=0.3):
    ax.fill_between(x, lo, hi, color=color, alpha=fill_alpha,
                    linestyle=linestyle)
    ax.plot(x, mid, color=color, linestyle=linestyle, label=label)


def plot_one(model: str, scotese_df):
    flux = cfg.load_carbon_flux(model)
    gast = cfg.load_global_temperature(model)
    times = flux.index.to_numpy()

    fig, ax, ax2 = _make_panels()

    # ------------------------------------------------------------------
    # Upper panel — net carbon outflux + total outflux, two rift variants
    # ------------------------------------------------------------------
    corrected_net   = flux["net_carbon_outflux_no_sed_unbiased_rifts"]
    uncorrected_net = flux["net_carbon_outflux_no_sed_biased_rifts"]
    corrected_tot   = flux["gross_atmospheric_influx_no_sed_unbiased_rifts"]
    uncorrected_tot = flux["gross_atmospheric_influx_no_sed_biased_rifts"]

    _fill_and_plot(ax, times,
                   corrected_net["mean"], corrected_net["min"], corrected_net["max"],
                   color="darkseagreen", label="Net carbon outflux (Corrected rift length)")
    _fill_and_plot(ax, times,
                   uncorrected_net["mean"], uncorrected_net["min"], uncorrected_net["max"],
                   color="darkseagreen", label="Net carbon outflux (Uncorrected rift length)",
                   linestyle=":", fill_alpha=0.2)
    _fill_and_plot(ax, times,
                   corrected_tot["mean"], corrected_tot["min"], corrected_tot["max"],
                   color="C4", label="Total outflux (Corrected rift length)")
    _fill_and_plot(ax, times,
                   uncorrected_tot["mean"], uncorrected_tot["min"], uncorrected_tot["max"],
                   color="C4", label="Total outflux (Uncorrected rift length)",
                   linestyle=":", fill_alpha=0.2)

    ax.plot(times, np.zeros_like(times), color="0.5")
    ax.set_xlim(170, 0)
    ax.set_xticks(np.arange(170, -1, -10), minor=True)
    ax.set_xticks(np.arange(170, -1, -50), minor=False)
    ax.set_yticks(np.arange(-50, 161, 5), minor=True)
    ax.tick_params(direction="in", which="major", length=5,
                   bottom=True, top=True, right=True)
    ax.tick_params(direction="in", which="minor", length=2.5,
                   bottom=True, top=True, right=True)
    ax.set_xticklabels([])  # bottom panel carries the labels
    ax.spines["bottom"].set_visible(True)
    ax.grid(alpha=0.1, which="both")
    ax.grid(alpha=0.3, which="major")
    ax.set_ylabel("Carbon\nflux (MtC/yr)", rotation="vertical", fontsize=12)
    ax.legend(fontsize=8, loc="upper left", ncol=2, frameon=False)

    # ------------------------------------------------------------------
    # Lower panel — GAST (Mills) with grey envelope + Scotese curve
    # ------------------------------------------------------------------
    g_t = gast.index.to_numpy()
    ax2.fill_between(g_t, gast["GAST_min"], gast["GAST_max"],
                     color="0.5", alpha=0.3, zorder=1)
    ax2.plot(g_t, gast["GAST_mean"], color="k", label="Mills", zorder=3)

    # Scotese — slice to the 170-0 Ma window so it doesn't dominate the y-axis
    scot = scotese_df[(scotese_df["age"] >= 0) & (scotese_df["age"] <= 170)]
    ax2.plot(scot["age"], scot["scotese_gast"], color="blue",
             label="Scotese", zorder=9)

    ax2.set_xlim(170, 0)
    ax2.set_xticks(np.arange(170, -1, -10), minor=True)
    ax2.set_xticks(np.arange(170, -1, -50), minor=False)
    ax2.minorticks_on()
    ax2.tick_params(direction="in", which="major", length=5,
                    bottom=True, top=True, right=True)
    ax2.tick_params(direction="in", which="minor", length=2.5,
                    bottom=True, top=True, right=True)
    ax2.set_xlabel("Age (Ma)", fontsize=12)
    ax2.set_ylabel("Global\ntemperature\n(℃)", rotation="vertical",
                   fontsize=12)
    ax2.legend(fontsize=9, loc="upper right", frameon=False)
    ax2.grid(alpha=0.1, which="both")

    fig.suptitle(f"170 -> 0 Ma: net & total atmospheric flux + GAST — {model}",
                 fontsize=11, y=1.02)
    cfg.save_figure(fig, f"170-0_diff_atmospheric_flux_gast_{model}")
    plt.close(fig)


def main():
    scotese_df = cfg.load_scotese_gast()
    for model in cfg.MODELS:
        print(f"[{model}]")
        plot_one(model, scotese_df)


if __name__ == "__main__":
    main()
