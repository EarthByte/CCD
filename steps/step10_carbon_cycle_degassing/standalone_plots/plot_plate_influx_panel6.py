#!/usr/bin/env python3
"""
Standalone reproduction of panel (f) of plate_influx_comparison.png from
Notebook 02 (02-Subducted-Carbon.ipynb), generated independently for the
DM2026 and BW1991 CCD scenarios.

Source: cell 46 of 02-Subducted-Carbon.ipynb (the c == 5 branch — all
reservoirs plotted together with the total in black).
Data:   02_plate_influx.csv (per model).

Outputs (under ./output/):
    plate_influx_panel6_DM2026.png/pdf
    plate_influx_panel6_BW1991.png/pdf
"""
from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import _common as cfg


# Component → CSV column → display label → colour (same as the notebook).
COMPONENTS = [
    ("sediments",          "Sediments",                "C0"),
    ("crust",              "Crust",                    "C2"),
    ("lithosphere",        "Lithosphere",              "C3"),
    # Serpentinite is split into two sub-series in the panel — bending (solid)
    # and MOR (dashed). Both share the same colour.
    ("serpentinite_bending", "Serpentinite (plate bending)", "C1"),
    ("serpentinite_mor",     "Serpentinite (MOR)",           "C1"),
    ("organic_sediments",  "Organic Sediment",         "C4"),
]


def plot_one(model: str):
    df = cfg.load_plate_influx(model)
    max_time = float(df.index.max())

    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    for col, label, colour in COMPONENTS:
        linestyle = "--" if col == "serpentinite_mor" else "-"
        ax.fill_between(df.index,
                        df[col]["min"].to_numpy(),
                        df[col]["max"].to_numpy(),
                        color=colour, alpha=0.2)
        ax.plot(df.index, df[col]["mean"].to_numpy(),
                color=colour, linestyle=linestyle, label=label)

    # Total influx — grey shading + black line.
    ax.fill_between(df.index,
                    df["total_influx"]["min"].to_numpy(),
                    df["total_influx"]["max"].to_numpy(),
                    color="0.5", alpha=0.5)
    ax.plot(df.index, df["total_influx"]["mean"].to_numpy(),
            color="k", label="Total")

    ax.set_xlim(max_time, 0)
    ax.set_xlabel("Age (Ma)")
    ax.set_ylabel("Rate of plate influx (Mt C/yr)")
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.tick_params(axis="both", direction="in", which="both",
                   right=True, top=True)
    ax.legend(frameon=False, loc="upper left", fontsize=9, ncols=1)
    ax.set_title(f"Plate influx — {model}", fontsize=11)

    cfg.save_figure(fig, f"plate_influx_panel6_{model}")
    plt.close(fig)


def main():
    for model in cfg.MODELS:
        print(f"[{model}]")
        plot_one(model)


if __name__ == "__main__":
    main()
