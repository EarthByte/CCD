#!/usr/bin/env python3
"""
Standalone reproduction of reservoir_flux.png from Notebook 05
(05-Atmospheric-Carbon.ipynb), one figure per CCD model.

Three series, all from 05_overriding_plate_storage.csv:
    Slab storage             — slab_storage_*           (crimson)
    Overriding plate storage — overriding_plate_storage_* (blue)
    Subduction outflux       — subduction_atmospheric_influx_* (C1)

Source: cell 41 of 05-Atmospheric-Carbon.ipynb.

Outputs (under ./output/):
    reservoir_flux_DM2026.png/pdf
    reservoir_flux_BW1991.png/pdf
"""
from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import _common as cfg


SERIES = [
    ("slab_storage",                    "Slab storage",             "crimson"),
    ("overriding_plate_storage",        "Overriding plate storage", "blue"),
    ("subduction_atmospheric_influx",   "Subduction outflux",       "C1"),
]


def plot_one(model: str):
    df = cfg.load_overriding_plate_storage(model)
    times = df.index.to_numpy()
    max_time = float(times.max())

    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    for stem, label, colour in SERIES:
        ax.fill_between(times,
                        df[f"{stem}_min"].to_numpy(),
                        df[f"{stem}_max"].to_numpy(),
                        color=colour, alpha=0.5)
        ax.plot(times, df[f"{stem}_mean"].to_numpy(),
                color=colour, label=label)

    ax.set_xlim(max_time, 0)
    ax.set_xlabel("Age (Ma)")
    ax.set_ylabel("Reservoir flux (Mt C/yr)")
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.set_ticks(np.arange(0., max_time + 1, 20.), minor=False)
    ax.tick_params(axis="both", direction="in", which="both",
                   right=True, top=True)
    ax.grid(alpha=0.1, which="minor")
    ax.grid(alpha=0.3, which="major")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title(f"Reservoir flux — {model}", fontsize=11)

    cfg.save_figure(fig, f"reservoir_flux_{model}")
    plt.close(fig)


def main():
    for model in cfg.MODELS:
        print(f"[{model}]")
        plot_one(model)


if __name__ == "__main__":
    main()
