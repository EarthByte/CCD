#!/usr/bin/env python3
"""
Net atmospheric carbon outflux, DM2026 (current run) — end-member carbonate
scenarios.

The growing pelagic carbonate reservoir bounds the net solid-Earth outflux
between two end-members:
  * carbonate reservoir is a genuinely NEW sink  -> count the full sediment
    reservoir on the storage side -> LOWER net outflux (too little degassing);
  * carbonate reservoir is only a shallow->deep REDISTRIBUTION of burial ->
    exclude the sediment reservoir from storage -> HIGHER net outflux (too much
    degassing).
The truth lies between the two. Five curves are shown:
  1. Gross atmospheric outflux                       (context)
  2. Upper-plate storage INCL. carbonate sediment    (sediments+serpentinite+crust)
  3. Upper-plate storage EXCL. carbonate sediment    (serpentinite+crust)
  4. Net outflux WITH sediment sink   (new-sink end-member, min degassing)
  5. Net outflux WITHOUT sediment sink (shift end-member, max degassing)

Sources (current Alfonso_etal_2024_DM26/Outputs run):
  gross / storage-incl / net-with  : 05_atmospheric_influx_all_sources.csv
  storage-excl                     : 02_plate_influx.csv (serpentinite_total + crust)
  net-without (net_carbon_outflux) : 05_net_carbon_outflux.csv

Output: ./output/net_atmospheric_influx_DM2026.png/pdf
"""
from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import _common as cfg


def _tight_ylim(values, pad_frac=0.05):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -1.0, 1.0
    lo, hi = float(finite.min()), float(finite.max())
    if hi == lo:
        return lo - 1.0, hi + 1.0
    pad = pad_frac * (hi - lo)
    lo = min(lo - pad, 0.0) if lo > -pad else lo - pad
    return lo, hi + pad


def plot_one(model: str):
    atm   = cfg.load_atmospheric_influx(model)
    plate = cfg.load_plate_influx(model).reindex(atm.index)
    net   = cfg.load_net_carbon_outflux(model).reindex(atm.index)
    t = atm.index.to_numpy(dtype=float)

    def flat(stem):
        return {q: atm[f"{stem}_{q}"].to_numpy() for q in ("min", "mean", "max")}

    gross       = flat("gross_atmospheric_outflux_unbiased_rift")
    stor_incl   = flat("gross_upper_plate_influx_with_sed")
    stor_excl   = {q: plate[("serpentinite_total", q)].to_numpy() + plate[("crust", q)].to_numpy()
                   for q in ("min", "mean", "max")}
    net_with    = flat("net_atmospheric_influx_unbiased_and_sed")               # new-sink end-member
    net_without = {q: net[("net_carbon_outflux", q)].to_numpy() for q in ("min", "mean", "max")}  # shift

    fig, ax = plt.subplots(figsize=(9, 6))
    pooled = []

    # envelope between the two net end-members (the plausible range)
    ax.fill_between(t, net_with["mean"], net_without["mean"], color="0.6", alpha=0.18, zorder=0,
                    label="Net-outflux range (end-members)")

    # 1. gross (context)
    ax.plot(t, gross["mean"], color="#8172B3", lw=1.6, label="Gross atmospheric outflux")
    ax.fill_between(t, gross["min"], gross["max"], color="#8172B3", alpha=0.12)
    # 2/3. storage incl / excl carbonate sediment
    ax.plot(t, stor_incl["mean"], color="#2A9D8F", lw=1.6, ls="-",
            label="Upper-plate storage (incl. carbonate sediment)")
    ax.plot(t, stor_excl["mean"], color="#2A9D8F", lw=1.6, ls="--",
            label="Upper-plate storage (excl. carbonate sediment)")
    # 4. net WITH sediment sink (new-sink end-member; min degassing)
    ax.plot(t, net_with["mean"], color="#4C72B0", lw=2.0,
            label="Net outflux — carbonate = new sink (min degassing)")
    ax.fill_between(t, net_with["min"], net_with["max"], color="#4C72B0", alpha=0.20)
    # 5. net WITHOUT sediment sink (shift end-member; max degassing)
    ax.plot(t, net_without["mean"], color="#C44E52", lw=2.0,
            label="Net outflux — carbonate = shallow→deep shift (max degassing)")
    ax.fill_between(t, net_without["min"], net_without["max"], color="#C44E52", alpha=0.20)

    for d in (gross, stor_incl, stor_excl, net_with, net_without):
        pooled.extend([d["min"], d["max"], d["mean"]])
    ax.plot(t, np.zeros_like(t), color="0.5", lw=0.8)

    ax.set_xlim(float(t.max()), float(t.min()))
    ax.set_ylim(*_tight_ylim(np.concatenate(pooled)))
    ax.set_xlabel("Time (Ma)")
    ax.set_ylabel("Carbon flux\n(Mt C/a)", rotation="horizontal")
    ax.yaxis.set_label_coords(-0.16, 0.45)
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.tick_params(direction="in", which="major", length=5, bottom=True, top=True, right=True)
    ax.tick_params(direction="in", which="minor", length=2.5, bottom=True, top=True, right=True)
    ax.grid(alpha=0.1, which="both"); ax.grid(alpha=0.3, which="major")
    ax.legend(loc="upper left", frameon=False, fontsize=7.6)
    ax.set_title(f"Net atmospheric carbon outflux — end-member carbonate scenarios — {model}",
                 fontsize=10)

    cfg.save_figure(fig, f"net_atmospheric_influx_{model}")
    plt.close(fig)


def main():
    for model in cfg.MODELS:
        print(f"[{model}]")
        plot_one(model)


if __name__ == "__main__":
    main()
