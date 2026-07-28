#!/usr/bin/env python3
"""
Standalone reproduction of panel (f) of carbon_subducted_comparison.png
from Notebook 02 — the panel with all carbon reservoirs combined and the
total in black. No 'f' panel label is drawn (per the user's request).

Source: cell 40 of 02-Subducted-Carbon.ipynb (the c == 5 branch).
Data:   02_subducted_carbon.csv (per model).

Outputs (under ./output/):
    carbon_subducted_panel6_DM2026.png/pdf
    carbon_subducted_panel6_BW1991.png/pdf
"""
from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import _common as cfg


# Component CSV column, label, colour. Matches Notebook 02 cell 40.
COMPONENTS = [
    ("sediments",            "Sediments",                  "C0"),
    ("crust",                "Crust",                      "C2"),
    ("lithosphere",          "Lithosphere",                "C3"),
    ("serpentinite_bending", "Serpentinite (slab bending)", "C1"),
    ("serpentinite_mor",     "Serpentinite (MOR)",          "C1"),
    ("organic_sediments",    "Organic Sediment",           "C4"),
]


def plot_one(model: str):
    df = cfg.load_subducted_carbon(model)
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

    # Total subducted — grey shading + black line.
    ax.fill_between(df.index,
                    df["total"]["min"].to_numpy(),
                    df["total"]["max"].to_numpy(),
                    color="0.5", alpha=0.5)
    ax.plot(df.index, df["total"]["mean"].to_numpy(),
            color="k", label="Total")

    ax.set_xlim(max_time, 0)
    ax.set_xlabel("Age (Ma)")
    ax.set_ylabel("Subducted carbon (Mt C/yr)")
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.tick_params(axis="both", direction="in", which="both",
                   right=True, top=True)
    ax.legend(frameon=False, loc="upper left", fontsize=9, ncols=1)
    ax.set_title(f"Subducted carbon — {model}", fontsize=11)

    # (Note: no 'f' panel label - explicitly omitted per the user's request.)
    cfg.save_figure(fig, f"carbon_subducted_panel6_{model}")
    plt.close(fig)


def main():
    for model in cfg.MODELS:
        print(f"[{model}]")
        plot_one(model)


if __name__ == "__main__":
    main()
