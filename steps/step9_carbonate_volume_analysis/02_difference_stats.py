#!/usr/bin/env python3
"""
=============================================================================
02_difference_stats.py  -  per-time summary stats + picked-times sidecar
=============================================================================

For each 1 Myr time slice in [0, 150] Ma compute frame-invariant summary
statistics of the DM2026 - BW1991 carbonate-thickness difference field over
the common-deposition mask (cells where BOTH products produce carbonate),
then pick four well-separated maxima of mean |delta| using greedy non-
maximum suppression with a 25 Myr exclusion radius.

Mirrors the pyBacktrack-paper Fig 10 selection so the writeup can use the
same statistical hook ("the four times where the products disagree most,
well-separated through the Cretaceous-Cenozoic record").

Outputs:
    output/stats/difference_stats.csv
        per-time: n_cells, mean, std, median, IQR, max_abs, mean_abs.
    output/stats/picked_times.txt
        the four greedy-NMS-picked ages and their mean_abs values.

Run:
    python 02_difference_stats.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as cfg


def per_time_stats(t: int) -> dict | None:
    """Single-row dict of stats at age t, or None if no overlap."""
    a = cfg.load_grid("DM2026", t).values
    b = cfg.load_grid("BW1991", t).values
    m = np.isfinite(a) & np.isfinite(b)
    if not m.any():
        return None
    d = (a - b)[m]
    abs_d = np.abs(d)
    q25, q75 = np.percentile(d, (25, 75))
    return dict(
        time_Ma   = t,
        n_cells   = int(m.sum()),
        mean      = float(d.mean()),
        std       = float(d.std()),
        median    = float(np.median(d)),
        iqr       = float(q75 - q25),
        max_abs   = float(abs_d.max()),
        mean_abs  = float(abs_d.mean()),
    )


def write_csv(rows, path: Path) -> None:
    cols = ("time_Ma", "n_cells",
            "mean", "std", "median", "iqr",
            "max_abs", "mean_abs")
    with open(path, "w") as fh:
        fh.write("# Per-time DM2026 - BW1991 carbonate-thickness statistics\n")
        fh.write("# over the cells where both products produce carbonate.\n")
        fh.write("# Positive = DM2026 produces more carbonate than BW1991.\n")
        fh.write(",".join(cols) + "\n")
        for r in rows:
            vals = []
            for c in cols:
                v = r[c]
                vals.append(str(v) if isinstance(v, int) else f"{v:.4f}")
            fh.write(",".join(vals) + "\n")
    print(f"  wrote {path}")


def pick_times(rows, n_picks: int, exclusion_myr: int) -> list[dict]:
    """Greedy non-maximum suppression on mean_abs over the time axis."""
    sorted_rows = sorted(rows, key=lambda r: r["mean_abs"], reverse=True)
    chosen: list[dict] = []
    for r in sorted_rows:
        if all(abs(r["time_Ma"] - c["time_Ma"]) >= exclusion_myr
               for c in chosen):
            chosen.append(r)
            if len(chosen) == n_picks:
                break
    return sorted(chosen, key=lambda r: r["time_Ma"])


def main() -> int:
    print(f"\nComputing per-time stats over {len(cfg.TIMES)} ages "
          f"({cfg.MIN_TIME_MA}..{cfg.MAX_TIME_MA} Ma) ...")
    rows = []
    for t in cfg.TIMES:
        s = per_time_stats(t)
        if s is None:
            print(f"  t={t:>3} Ma  -- no common-mask overlap, skipping")
            continue
        rows.append(s)
        if t % 10 == 0:
            print(f"  t={t:>3} Ma | "
                  f"n={s['n_cells']:>7d} | "
                  f"mean_abs={s['mean_abs']:6.2f} m | "
                  f"max_abs={s['max_abs']:7.1f} m")

    write_csv(rows, cfg.OUTPUT_DIR / "stats" / "difference_stats.csv")

    picks = pick_times(rows, cfg.N_PICKED_TIMES, cfg.NMS_EXCLUSION_MYR)
    print(f"\nPicked {len(picks)} well-separated times "
          f"(greedy NMS, >= {cfg.NMS_EXCLUSION_MYR} Myr apart):")
    for p in picks:
        print(f"  {p['time_Ma']:>3} Ma | "
              f"mean_abs={p['mean_abs']:.2f} m, "
              f"max_abs={p['max_abs']:.1f} m, "
              f"iqr={p['iqr']:.2f} m")

    out_pick = cfg.OUTPUT_DIR / "stats" / "picked_times.txt"
    with open(out_pick, "w") as fh:
        fh.write("# Picked times for the maximum-disagreement analysis.\n")
        fh.write("# Method: greedy NMS on per-time mean(|DM2026 - BW1991|),\n")
        fh.write(f"# exclusion radius = {cfg.NMS_EXCLUSION_MYR} Myr, "
                 f"n_picks = {cfg.N_PICKED_TIMES}.\n")
        fh.write("time_Ma,mean_abs_m,max_abs_m,iqr_m\n")
        for p in picks:
            fh.write(f"{p['time_Ma']},{p['mean_abs']:.4f},"
                     f"{p['max_abs']:.4f},{p['iqr']:.4f}\n")
    print(f"  wrote {out_pick}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
