#!/usr/bin/env python3
"""
=============================================================================
04_region_characterisation.py
=============================================================================

Where, in geographic terms, do the two CCD curves disagree most through
time?

For every 1 Myr slice in [0, 150] Ma we bin the |DM2026 - BW1991| field
into a (latitude-band x ocean-basin) matrix and compute the mean per
matrix cell.  Latitude bands are 10 deg wide (-90..+90 -> 18 bands);
basins are present-day Pacific / Atlantic / Indian / Arctic / Southern
(see ``config.basin_of``).  This gives a (n_times x 18*5 = 90) wide
matrix that we serialise as CSV and visualise as a heatmap.

Three deliverables:

    output/stats/regional_mean_abs.csv
        wide CSV: rows = time, cols = (basin, lat_band_centre).

    output/figures/04_region_heatmap.png / pdf
        heatmap of mean |delta| per (basin x lat-band) over time.

    output/stats/picked_time_region_tables.csv
        for each of the four picked times (read from
        output/stats/picked_times.txt) the top-3 (basin, lat band) cells
        ranked by mean |delta|.  Goes into the writeup.

Run AFTER 02_difference_stats.py (so the picked-times sidecar exists).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as cfg

# Crameri's batlow is the publication cmap; fall back to viridis (also
# perceptually uniform, matplotlib built-in) if cmcrameri isn't installed.
try:
    import cmcrameri.cm as _cmc           # noqa: F401  (registers cmaps)
    HEATMAP_CMAP = "cmc.batlow"
except ImportError:
    HEATMAP_CMAP = "viridis"
    print("  [info] cmcrameri not installed; using matplotlib's viridis "
          "as the heatmap cmap. `pip install cmcrameri` for the Crameri "
          "batlow used in the publication figures.")


PICKED_TIMES_PATH = cfg.OUTPUT_DIR / "stats" / "picked_times.txt"

TOP_N_PER_PICK = 3       # top-N (basin, lat-band) cells reported per picked time


# ----------------------------------------------------------------------------
# Build the per-cell (basin, lat-band-index) classifier once -- it's the
# same for every time slice because the rasters share the present-day grid.
# ----------------------------------------------------------------------------
def build_classifier(lat_vec: np.ndarray, lon_vec: np.ndarray):
    """Return (basin_2d, lat_band_idx_2d) — both shape (n_lat, n_lon).

    `basin_2d` holds basin strings, `lat_band_idx_2d` holds the index into
    cfg.LAT_BAND_CENTRES (-1 for cells outside any band).
    """
    lat_2d, lon_2d = np.meshgrid(lat_vec, lon_vec, indexing="ij")
    basin_2d = cfg.basin_of(lat_2d, lon_2d)
    band_idx = np.digitize(lat_2d, cfg.LAT_BAND_EDGES) - 1
    band_idx = np.where(
        (band_idx >= 0) & (band_idx < len(cfg.LAT_BAND_CENTRES)),
        band_idx, -1)
    return basin_2d, band_idx


def regional_aggregate(t: int, basin_2d, band_idx) -> np.ndarray:
    """Return a (n_basins, n_lat_bands) matrix of mean |delta| at age t.

    NaN for empty cells.
    """
    a = cfg.load_grid("DM2026", t).values
    b = cfg.load_grid("BW1991", t).values
    abs_d = np.abs(a - b)
    valid = np.isfinite(a) & np.isfinite(b) & (band_idx >= 0)

    nb, nl = len(cfg.BASINS), len(cfg.LAT_BAND_CENTRES)
    out = np.full((nb, nl), np.nan, dtype=float)
    for bi, name in enumerate(cfg.BASINS):
        mask_b = valid & (basin_2d == name)
        if not mask_b.any():
            continue
        rows = abs_d[mask_b]
        bands = band_idx[mask_b]
        for li in range(nl):
            sel = bands == li
            if sel.any():
                out[bi, li] = float(rows[sel].mean())
    return out


def load_picked_times() -> list[int]:
    if not PICKED_TIMES_PATH.exists():
        print(f"  [warn] {PICKED_TIMES_PATH} missing -- "
              "run 02_difference_stats.py first")
        return []
    picks = []
    with open(PICKED_TIMES_PATH) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("time_Ma"):
                continue
            picks.append(int(line.split(",")[0]))
    return picks


def main() -> int:
    picked_times = load_picked_times()
    print(f"  picked-times: {picked_times}")

    ref = cfg.load_grid("DM2026", 0)
    lat_vec = ref["lat"].values
    lon_vec = ref["lon"].values
    basin_2d, band_idx = build_classifier(lat_vec, lon_vec)

    # ---- per-time aggregation
    print(f"\nAggregating per-time (basin x lat-band) means over "
          f"{len(cfg.TIMES)} ages ...")
    mats = []                                        # list of (nb, nl) arrays
    for t in cfg.TIMES:
        mats.append(regional_aggregate(t, basin_2d, band_idx))
        if t % 10 == 0:
            print(f"  t={t:>3} Ma  done")

    nb, nl = len(cfg.BASINS), len(cfg.LAT_BAND_CENTRES)
    cube = np.stack(mats, axis=0)                    # (n_times, nb, nl)

    # ---- wide CSV: rows=time, cols=(basin, band_centre)
    wide_path = cfg.OUTPUT_DIR / "stats" / "regional_mean_abs.csv"
    header = ["time_Ma"]
    for bn in cfg.BASINS:
        for c in cfg.LAT_BAND_CENTRES:
            header.append(f"{bn}_{int(c):+d}")
    with open(wide_path, "w") as fh:
        fh.write("# Mean |DM2026 - BW1991| (m) per (basin, 10-deg lat band).\n")
        fh.write("# Lat-band tag = band centre in degrees (signed).\n")
        fh.write(",".join(header) + "\n")
        for ti, t in enumerate(cfg.TIMES):
            row = [str(t)]
            for bi in range(nb):
                for li in range(nl):
                    v = cube[ti, bi, li]
                    row.append(f"{v:.4f}" if np.isfinite(v) else "")
            fh.write(",".join(row) + "\n")
    print(f"  wrote {wide_path}")

    # ---- heatmap figure
    # Reshape to (nb*nl) rows x n_times cols, then plot with rows
    # grouped by basin and labelled by band centre.
    H = cube.transpose(1, 2, 0).reshape(nb * nl, len(cfg.TIMES))  # (nb*nl, n_times)

    fig, ax = plt.subplots(figsize=(13, 10))
    # Mask NaNs as transparent
    masked = np.ma.masked_invalid(H)
    vmax = float(np.nanpercentile(H, 99))
    vmin = max(0.1, float(np.nanpercentile(H[H > 0], 1)) if (H > 0).any() else 0.1)
    im = ax.imshow(masked, aspect="auto", origin="lower",
                   cmap=HEATMAP_CMAP,
                   norm=LogNorm(vmin=vmin, vmax=vmax),
                   extent=[cfg.TIMES[0] - 0.5, cfg.TIMES[-1] + 0.5,
                           -0.5, nb * nl - 0.5])
    ax.set_xlim(cfg.MAX_TIME_MA, cfg.MIN_TIME_MA)
    ax.set_xlabel("Age (Ma)", fontsize=13)
    ax.set_ylabel("Region (basin / latitude-band centre)", fontsize=13)

    # Row labels grouped by basin
    yticks_minor = []
    yticks_minor_labels = []
    yticks_major = []
    yticks_major_labels = []
    for bi, bn in enumerate(cfg.BASINS):
        for li in range(nl):
            yticks_minor.append(bi * nl + li)
            yticks_minor_labels.append(f"{int(cfg.LAT_BAND_CENTRES[li]):+d}")
        yticks_major.append(bi * nl + nl / 2 - 0.5)
        yticks_major_labels.append(bn)
        ax.axhline(bi * nl - 0.5, color="white", lw=0.5)
    ax.set_yticks(yticks_minor, minor=True)
    ax.set_yticklabels(yticks_minor_labels, minor=True, fontsize=7)
    ax.set_yticks(yticks_major)
    ax.set_yticklabels(yticks_major_labels, fontsize=11, fontweight="bold")
    ax.tick_params(axis="y", which="major", length=0, pad=40)
    ax.tick_params(axis="y", which="minor", length=4)

    # Highlight picked times with vertical lines
    for t in picked_times:
        ax.axvline(t, color="red", lw=1.0, linestyle="--", alpha=0.8)
        ax.text(t, nb * nl - 0.5, f"{t} Ma", color="red",
                ha="center", va="bottom", fontsize=10, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("mean |DM2026 - BW1991| (m, log scale)", fontsize=11)
    ax.set_title("Where the two CCD curves disagree most -- "
                 "regional time-series of mean |delta thickness|",
                 fontsize=13)
    fig.tight_layout()
    out_png = cfg.OUTPUT_DIR / "figures" / "04_region_heatmap.png"
    out_pdf = cfg.OUTPUT_DIR / "figures" / "04_region_heatmap.pdf"
    fig.savefig(out_png, dpi=200)
    fig.savefig(out_pdf)
    print(f"  wrote {out_png}")
    print(f"  wrote {out_pdf}")
    plt.close(fig)

    # ---- per-picked-time top-N tables
    if picked_times:
        table_path = cfg.OUTPUT_DIR / "stats" / "picked_time_region_tables.csv"
        with open(table_path, "w") as fh:
            fh.write("# Top-{n} (basin, lat-band) cells ranked by mean "
                     "|DM2026 - BW1991| at each picked time.\n".format(
                         n=TOP_N_PER_PICK))
            fh.write("time_Ma,rank,basin,lat_band_centre_deg,mean_abs_m\n")
            for t in picked_times:
                ti = cfg.TIMES.index(t)
                flat = []
                for bi in range(nb):
                    for li in range(nl):
                        v = cube[ti, bi, li]
                        if np.isfinite(v):
                            flat.append((v, cfg.BASINS[bi],
                                         int(cfg.LAT_BAND_CENTRES[li])))
                flat.sort(reverse=True)
                for rank, (v, bn, c) in enumerate(flat[:TOP_N_PER_PICK], 1):
                    fh.write(f"{t},{rank},{bn},{c:+d},{v:.4f}\n")
        print(f"  wrote {table_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
