#!/usr/bin/env python3
"""
Run every plotting script in this folder. Equivalent to executing each
plot_*.py individually but cheaper because pandas/matplotlib are imported
just once.
"""
from __future__ import annotations

import importlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SCRIPTS = [
    "plot_plate_influx_panel6",
    "plot_carbon_subducted_panel6",
    "plot_reservoir_flux",
    "plot_net_atmospheric_influx",
    "plot_atmosphere_influx_comparison",
    "plot_diff_atmospheric_flux_gast",
    "plot_diff_atmospheric_flux_gast_ratschbacher",
]


def main():
    for name in SCRIPTS:
        print(f"\n===== {name} =====")
        mod = importlib.import_module(name)
        mod.main()


if __name__ == "__main__":
    main()
