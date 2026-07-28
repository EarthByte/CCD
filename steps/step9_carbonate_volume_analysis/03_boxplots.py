#!/usr/bin/env python3
"""
=============================================================================
03_boxplots.py  -  2-panel comparison figure
=============================================================================

Panel (a): xy time-series of compacted carbonate thickness over 0-150 Ma.
           Plots DM2026 and BW1991 means with the +/- 1 standard-deviation
           envelope shaded around each line.  Statistics are computed per
           1 Myr slice over the common deposition mask.

Panel (b): signed delta = DM2026 - BW1991 distribution in 5 Myr bins.
           Same boxplot construction as Fig 11 of the pyBacktrack repo.
           The four picked times of largest mean |delta| (from
           02_difference_stats.py) are highlighted in red.

This replaces the earlier 3-panel layout where the DM2026 and BW1991
boxplot panels were too visually similar to compare at a glance.

Outputs:
    output/figures/03_thickness_and_delta.png
    output/figures/03_thickness_and_delta.pdf
    output/stats/thickness_mean_std.csv      (per-1Myr means + std for both)
    output/stats/boxplot_delta_summary.csv   (per-5Myr-bin delta stats)
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as cfg


# ----------------------------------------------------------------------------
# Tunables
# ----------------------------------------------------------------------------
SUBSAMPLE_MAX_CELLS = 200_000      # cap per-bin boxplot sample size
RNG = np.random.default_rng(20260603)

AXIS_LABEL_SIZE = 14
TICK_LABEL_SIZE = 12
TITLE_SIZE = 14
LEGEND_SIZE = 12

DM_COLOR = "#1f77b4"            # blue
BW_COLOR = "#d62728"            # red

PICKED_TIMES_PATH = cfg.OUTPUT_DIR / "stats" / "picked_times.txt"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def load_picked_times() -> list[int]:
    if not PICKED_TIMES_PATH.exists():
        print(f"  [warn] {PICKED_TIMES_PATH} missing -- run "
              "02_difference_stats.py first; no boxes will be highlighted")
        return []
    picks = []
    with open(PICKED_TIMES_PATH) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("time_Ma"):
                continue
            picks.append(int(line.split(",")[0]))
    return picks


def per_time_thickness_stats():
    """Per-1-Myr mean + std of compacted thickness over the common mask,
    for each source.  Returns three aligned 1-D arrays:

        times (n,)        ages in Ma
        dm    (n, 2)      cols (mean, std) for DM2026
        bw    (n, 2)      cols (mean, std) for BW1991
    """
    rows_dm = []
    rows_bw = []
    times = []
    print(f"\nPer-time thickness stats over {len(cfg.TIMES)} ages ...")
    for t in cfg.TIMES:
        a = cfg.load_grid("DM2026", t).values
        b = cfg.load_grid("BW1991", t).values
        m = np.isfinite(a) & np.isfinite(b)
        if not m.any():
            continue
        av = a[m]; bv = b[m]
        times.append(t)
        rows_dm.append((float(av.mean()), float(av.std())))
        rows_bw.append((float(bv.mean()), float(bv.std())))
        if t % 10 == 0:
            print(f"  t={t:>3} Ma | DM mean={rows_dm[-1][0]:6.1f}+/-{rows_dm[-1][1]:5.1f}"
                  f" | BW mean={rows_bw[-1][0]:6.1f}+/-{rows_bw[-1][1]:5.1f}")
    return (np.array(times),
            np.array(rows_dm, dtype=float),
            np.array(rows_bw, dtype=float))


def collect_delta_bin(t_lo: int, t_hi: int) -> np.ndarray:
    """Pool every cell delta value in [t_lo, t_hi) Ma into a 1-D array."""
    chunks = []
    for t in range(t_lo, t_hi):
        if t < cfg.MIN_TIME_MA or t > cfg.MAX_TIME_MA:
            continue
        a = cfg.load_grid("DM2026", t).values
        b = cfg.load_grid("BW1991", t).values
        m = np.isfinite(a) & np.isfinite(b)
        if m.any():
            chunks.append((a - b)[m])
    if not chunks:
        return np.array([])
    v = np.concatenate(chunks)
    if v.size > SUBSAMPLE_MAX_CELLS:
        v = v[RNG.choice(v.size, SUBSAMPLE_MAX_CELLS, replace=False)]
    return v


def summarise(v: np.ndarray) -> dict:
    if v.size == 0:
        return dict(n=0, mean=np.nan, median=np.nan, std=np.nan, iqr=np.nan)
    q25, q75 = np.percentile(v, (25, 75))
    return dict(n=int(v.size),
                mean=float(v.mean()),
                median=float(np.median(v)),
                std=float(v.std()),
                iqr=float(q75 - q25))


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> int:
    picked_times = load_picked_times()
    print(f"  picked-times for highlight: {picked_times}")

    # ---- Panel (a) data: per-1Myr mean +/- std for both products
    times, dm, bw = per_time_thickness_stats()

    # ---- write thickness_mean_std.csv
    thk_csv = cfg.OUTPUT_DIR / "stats" / "thickness_mean_std.csv"
    with open(thk_csv, "w") as fh:
        fh.write("# Per-1Myr mean +/- std of compacted carbonate thickness (m)\n")
        fh.write("# computed over the common deposition mask.\n")
        fh.write("time_Ma,mean_DM,std_DM,mean_BW,std_BW\n")
        for i, t in enumerate(times):
            fh.write(f"{int(t)},{dm[i,0]:.4f},{dm[i,1]:.4f},"
                     f"{bw[i,0]:.4f},{bw[i,1]:.4f}\n")
    print(f"  wrote {thk_csv}")

    # ---- Panel (b) data: 5 Myr-binned pooled delta values + summary
    edges = cfg.BIN_EDGES
    centres = cfg.BIN_CENTRES
    print(f"\nPooling delta cells per {cfg.BIN_WIDTH_MYR} Myr bin ...")
    delta_pools = []
    summary_rows = []
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        v = collect_delta_bin(int(lo), int(hi))
        delta_pools.append(v)
        s = summarise(v)
        summary_rows.append(dict(bin_centre_Ma=float(centres[i]), **s))
        print(f"  bin {int(lo):>3}-{int(hi):>3} Ma  n={s['n']:>7d}  "
              f"mean={s['mean']:+6.2f}  med={s['median']:+6.2f}  "
              f"iqr={s['iqr']:5.2f}")

    box_csv = cfg.OUTPUT_DIR / "stats" / "boxplot_delta_summary.csv"
    with open(box_csv, "w") as fh:
        fh.write("# 5 Myr-binned signed-delta (DM2026 - BW1991) stats.\n")
        fh.write("bin_centre_Ma,n,mean,median,std,iqr\n")
        for r in summary_rows:
            fh.write(f"{r['bin_centre_Ma']},{r['n']},"
                     f"{r['mean']:.4f},{r['median']:.4f},"
                     f"{r['std']:.4f},{r['iqr']:.4f}\n")
    print(f"  wrote {box_csv}")

    # ---- Figure
    fig, (ax_ts, ax_bx) = plt.subplots(
        nrows=2, ncols=1, figsize=(11, 9), sharex=True,
        gridspec_kw=dict(height_ratios=[1.0, 1.0]),
    )

    # Panel (a): mean +/- 1 sigma envelope, both sources overlaid
    for name, color, stats in (("DM2026", DM_COLOR, dm),
                               ("BW1991", BW_COLOR, bw)):
        mu = stats[:, 0]
        sd = stats[:, 1]
        ax_ts.fill_between(times, mu - sd, mu + sd,
                           color=color, alpha=0.20, linewidth=0,
                           label=f"{name} (mean +/- 1 sigma)")
        ax_ts.plot(times, mu, color=color, lw=2.0,
                   label=f"{name} mean")

    ax_ts.set_xlim(cfg.MAX_TIME_MA, cfg.MIN_TIME_MA)
    ax_ts.set_autoscalex_on(False)
    ax_ts.set_ylim(bottom=0)
    ax_ts.set_xticks(np.arange(cfg.MIN_TIME_MA, cfg.MAX_TIME_MA + 1, 10))
    ax_ts.set_xticks(np.arange(cfg.MIN_TIME_MA, cfg.MAX_TIME_MA + 1, 5),
                     minor=True)
    ax_ts.set_xticklabels([str(x) for x in
                           np.arange(cfg.MIN_TIME_MA, cfg.MAX_TIME_MA + 1, 10)])
    ax_ts.tick_params(axis="both", which="major", labelsize=TICK_LABEL_SIZE)
    ax_ts.set_ylabel("compacted thickness (m)", fontsize=AXIS_LABEL_SIZE)
    ax_ts.grid(True, axis="y", which="major", alpha=0.3, linestyle=":")
    ax_ts.legend(loc="upper right", fontsize=LEGEND_SIZE, framealpha=0.9, ncol=2)
    ax_ts.set_title("(a) Compacted carbonate sediment thickness "
                    "(mean +/- 1 sigma, common deposition mask)",
                    fontsize=TITLE_SIZE, loc="left")
    # mark picked times for cross-reference
    for t in picked_times:
        ax_ts.axvline(t, color="black", lw=0.6, linestyle=":", alpha=0.4)

    # Panel (b): delta boxplot, 5 Myr bins
    positions = centres
    width = cfg.BIN_WIDTH_MYR * 0.7
    bp = ax_bx.boxplot(
        delta_pools,
        positions=positions,
        widths=width,
        showfliers=False,
        patch_artist=True,
        medianprops=dict(color="black", linewidth=1.5),
        whiskerprops=dict(linewidth=0.7),
        capprops=dict(linewidth=0.7),
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("#fdbf6f")
        patch.set_edgecolor("black")
        patch.set_linewidth(0.5)
    if picked_times:
        for t in picked_times:
            i = int(np.argmin(np.abs(centres - t)))
            bp["boxes"][i].set_facecolor("#e31a1c")
            bp["boxes"][i].set_edgecolor("black")

    ax_bx.axhline(0, color="k", lw=0.5, linestyle="--")
    ax_bx.set_xlim(cfg.MAX_TIME_MA, cfg.MIN_TIME_MA)
    ax_bx.set_autoscalex_on(False)
    ax_bx.set_xticks(np.arange(cfg.MIN_TIME_MA, cfg.MAX_TIME_MA + 1, 10))
    ax_bx.set_xticks(np.arange(cfg.MIN_TIME_MA, cfg.MAX_TIME_MA + 1, 5),
                     minor=True)
    ax_bx.set_xticklabels([str(x) for x in
                           np.arange(cfg.MIN_TIME_MA, cfg.MAX_TIME_MA + 1, 10)])
    ax_bx.tick_params(axis="both", which="major", labelsize=TICK_LABEL_SIZE)
    ax_bx.set_xlabel("Age (Ma)", fontsize=AXIS_LABEL_SIZE)
    ax_bx.set_ylabel("thickness anomaly (m)", fontsize=AXIS_LABEL_SIZE)
    ax_bx.grid(True, axis="y", which="major", alpha=0.3, linestyle=":")
    ax_bx.set_title("(b) Signed difference DM2026 - BW1991 (5 Myr bins)",
                    fontsize=TITLE_SIZE, loc="left")

    # Sign-convention labels
    sign_x = 30
    trans = blended_transform_factory(ax_bx.transData, ax_bx.transAxes)
    ax_bx.text(sign_x, 0.97, "DM2026 thicker", transform=trans,
               ha="center", va="top", fontsize=LEGEND_SIZE,
               bbox=dict(facecolor="white", edgecolor="black",
                         boxstyle="round,pad=0.3", linewidth=0.5))
    ax_bx.text(sign_x, 0.03, "BW1991 thicker", transform=trans,
               ha="center", va="bottom", fontsize=LEGEND_SIZE,
               bbox=dict(facecolor="white", edgecolor="black",
                         boxstyle="round,pad=0.3", linewidth=0.5))

    fig.tight_layout()
    out_png = cfg.OUTPUT_DIR / "figures" / "03_thickness_and_delta.png"
    out_pdf = cfg.OUTPUT_DIR / "figures" / "03_thickness_and_delta.pdf"
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    print(f"  wrote {out_png}")
    print(f"  wrote {out_pdf}")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
