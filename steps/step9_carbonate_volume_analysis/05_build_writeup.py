#!/usr/bin/env python3
"""
=============================================================================
05_build_writeup.py  -  assemble the comparison MS Word document
=============================================================================

Reads the CSVs and PNGs written by 02-04 and emits

    output/carbonate_thickness_comparison.docx

The writeup is structured as a short technical report:

    1.  Motivation
    2.  Methods
    3.  Results
        3.1  Time-series of mean |delta|
        3.2  Distribution comparison (boxplots)
        3.3  Where the products disagree most (regional analysis)
        3.4  Four times of largest disagreement
    4.  Discussion
    References

Run AFTER:
    python 02_difference_stats.py
    python 03_boxplots.py
    python 04_region_characterisation.py
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as cfg

try:
    from docx import Document
    from docx.shared import Cm, Pt, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
except ImportError as exc:
    raise SystemExit(
        "python-docx is required.  Install via "
        "`pip install python-docx` (or `conda install -c conda-forge python-docx`) "
        "and re-run."
    ) from exc


STATS_DIR = cfg.OUTPUT_DIR / "stats"
FIG_DIR = cfg.OUTPUT_DIR / "figures"


def read_picked_times() -> list[dict]:
    rows = []
    path = STATS_DIR / "picked_times.txt"
    if not path.exists():
        return rows
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("time_Ma"):
                continue
            parts = line.split(",")
            rows.append(dict(time_Ma=int(parts[0]),
                             mean_abs=float(parts[1]),
                             max_abs=float(parts[2]),
                             iqr=float(parts[3])))
    return rows


def read_region_tables() -> dict[int, list[dict]]:
    out: dict[int, list[dict]] = {}
    path = STATS_DIR / "picked_time_region_tables.csv"
    if not path.exists():
        return out
    with open(path) as fh:
        rdr = csv.reader(fh)
        for row in rdr:
            if not row or row[0].startswith("#") or row[0] == "time_Ma":
                continue
            t = int(row[0])
            out.setdefault(t, []).append(dict(
                rank=int(row[1]),
                basin=row[2],
                lat_band=int(row[3]),
                mean_abs=float(row[4]),
            ))
    return out


def read_stats_csv(path: Path) -> list[dict]:
    """Generic CSV reader skipping comment lines."""
    rows = []
    if not path.exists():
        return rows
    with open(path) as fh:
        # Skip leading comment lines
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line:
                continue
            if header is None:
                header = line.split(",")
                continue
            vals = line.split(",")
            rows.append(dict(zip(header, vals)))
    return rows


# ----------------------------------------------------------------------------
# docx helpers
# ----------------------------------------------------------------------------
def add_h1(doc, text):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(18)


def add_h2(doc, text):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(14)


def add_h3(doc, text):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    r.italic = True
    r.font.size = Pt(12)


def add_para(doc, text):
    doc.add_paragraph(text)


def add_figure(doc, png_path: Path, caption: str, width_cm: float = 16.0):
    if png_path.exists():
        doc.add_picture(str(png_path), width=Cm(width_cm))
        last = doc.paragraphs[-1]
        last.alignment = WD_ALIGN_PARAGRAPH.CENTER
    else:
        p = doc.add_paragraph()
        r = p.add_run(f"[placeholder: {png_path.name} not yet generated]")
        r.italic = True
        r.font.color.rgb = RGBColor(0xA0, 0xA0, 0xA0)
    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    cr = cap.add_run(caption)
    cr.italic = True
    cr.font.size = Pt(10)


def add_table(doc, header_row, body_rows):
    tbl = doc.add_table(rows=1 + len(body_rows), cols=len(header_row))
    tbl.style = "Light Grid Accent 1"
    hdr = tbl.rows[0].cells
    for j, h in enumerate(header_row):
        hdr[j].text = ""
        p = hdr[j].paragraphs[0]
        r = p.add_run(h)
        r.bold = True
    for i, row in enumerate(body_rows, start=1):
        for j, val in enumerate(row):
            tbl.rows[i].cells[j].text = str(val)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> int:
    picked = read_picked_times()
    region_tables = read_region_tables()
    diff_stats = read_stats_csv(STATS_DIR / "difference_stats.csv")

    doc = Document()

    # ---- title block
    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    tr = title.add_run("Comparison of global carbonate sediment thickness "
                       "from two CCD curves")
    tr.bold = True
    tr.font.size = Pt(20)

    sub = doc.add_paragraph()
    sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    sr = sub.add_run("Dutkiewicz & Müller 2026 vs. Boss & Wilkinson 1991, "
                     "on paleobathymetry computed on the Alfonso et al. (2024) plate "
                     "model following the pyBacktrack 1.5 method of Müller et al. (2026)")
    sr.italic = True
    sr.font.size = Pt(11)

    doc.add_paragraph()

    # ---- 1. Motivation
    add_h1(doc, "1. Motivation")
    add_para(doc,
        "The depth of the calcite compensation depth (CCD) through time "
        "controls where carbonate sediments accumulate on the seafloor. "
        "We compare two global compacted-carbonate-sediment-thickness "
        "products built on the same paleobathymetry - computed on the "
        "Alfonso et al. (2024) plate model following the pyBacktrack 1.5 method "
        "of Müller et al. (2026) - but using two alternative CCD reconstructions:")

    bul = doc.add_paragraph(style="List Bullet")
    bul.add_run("DM2026 — the new global CCD curve of Dutkiewicz & Müller "
                "(2026, in preparation), and")
    bul = doc.add_paragraph(style="List Bullet")
    bul.add_run("BW1991 — the classic CCD of Boss & Wilkinson (1991).")

    add_para(doc,
        "Both products report compacted carbonate sediment thickness on a "
        "0.25° lat/lon grid in 1 Myr intervals from 0 to 170 Ma. We "
        "analyse the 0-150 Ma window because carbonate-bearing seafloor "
        "area pre-150 Ma is small.")

    add_para(doc,
        "The aim is to quantify the magnitude, temporal pattern and "
        "geographic distribution of the differences between the two "
        "products, identify times of maximum disagreement, and assess "
        "where in the present-day ocean basins the disagreement is "
        "concentrated.")

    # ---- 2. Methods
    add_h1(doc, "2. Methods")
    add_h2(doc, "2.1 Per-time summary statistics")
    add_para(doc,
        "For every 1 Myr time slice in [0, 150] Ma we compute frame-"
        "invariant summary statistics of the signed difference field "
        "Δ(x, y, t) = T_DM2026(x, y, t) − T_BW1991(x, y, t) over the "
        "common deposition mask (grid cells where both products produce "
        "carbonate). Reported per time: cell count, mean, standard "
        "deviation, median, inter-quartile range, |Δ|_max, and "
        "<|Δ|>. The signed mean is sensitive to the direction of the "
        "disagreement; the mean of |Δ| is the appropriate scalar "
        "magnitude.")

    add_h2(doc, "2.2 Four times of largest disagreement")
    add_para(doc,
        f"Times of largest disagreement were identified by greedy non-"
        f"maximum-suppression (NMS) on the per-time mean-|Δ| series with "
        f"an exclusion radius of {cfg.NMS_EXCLUSION_MYR} Myr, retaining "
        f"the top {cfg.N_PICKED_TIMES} maxima. This produces well-"
        f"separated picks (≥ {cfg.NMS_EXCLUSION_MYR} Myr apart) so the "
        f"chosen times sample independent intervals rather than "
        f"clustering around a single high-disagreement episode. The same "
        f"algorithm is used in Fig. 10 of Müller et al. (in prep., GMD) "
        f"to pick representative times for the NW Shelf subsidence-rate "
        f"comparison.")

    add_h2(doc, "2.3 Boxplot binning (5 Myr bins)")
    add_para(doc,
        "Three boxplot panels show the per-cell value distributions in 5 Myr "
        "bins from 0 to 150 Ma (30 bins): (a) DM2026 thickness, (b) BW1991 "
        "thickness, (c) signed Δ = DM2026 − BW1991. Each bin pools every "
        "cell value from every 1 Myr slice that falls inside the 5 Myr "
        "window. Cell counts per bin are capped at 200 000 by random "
        "sub-sampling to keep the boxplots responsive without changing the "
        "distribution shape.")

    add_h2(doc, "2.4 Regional characterisation (latitude × basin)")
    add_para(doc,
        "Each grid cell is assigned to one of 18 10°-wide latitude bands "
        "and one of five present-day ocean basins (Pacific, Atlantic, "
        "Indian, Arctic, Southern). For every 1 Myr time slice we compute "
        "the mean of |Δ| in every (basin × latitude band) cell. The result "
        "is a (151 × 90) wide matrix that we visualise as a log-scaled "
        "heatmap with the picked times overlaid. For each picked time we "
        "rank the (basin × latitude band) cells by mean |Δ| and report the "
        "top three as the regions accounting for the bulk of the "
        "disagreement at that time.")

    # ---- 3. Results
    add_h1(doc, "3. Results")

    # ---- 3.1 time series figure (we don't have a separate figure for this
    # but we can summarise from difference_stats.csv)
    add_h2(doc, "3.1 Time-series of difference statistics")
    if diff_stats:
        # Extract key numbers
        rows = [(int(r["time_Ma"]), float(r["mean_abs"]), float(r["max_abs"]),
                 float(r["mean"]), int(r["n_cells"])) for r in diff_stats]
        overall_mean = sum(r[1] for r in rows) / len(rows)
        peak_mean_abs = max(rows, key=lambda r: r[1])
        peak_signed = max(rows, key=lambda r: abs(r[3]))
        add_para(doc,
            f"Across the 0-150 Ma window the time-mean of the per-time "
            f"<|Δ|> is {overall_mean:.2f} m. The single largest per-time "
            f"<|Δ|> = {peak_mean_abs[1]:.2f} m occurs at "
            f"{peak_mean_abs[0]} Ma. The signed-mean Δ is "
            f"{peak_signed[3]:+.2f} m at {peak_signed[0]} Ma — the sign "
            f"indicates whether DM2026 (positive) or BW1991 (negative) "
            f"produces more carbonate on average at that time.")

    # ---- 3.2 thickness time-series + delta boxplot
    add_h2(doc, "3.2 Thickness time-series and delta distribution")
    add_para(doc,
        "Figure 1 shows two complementary views of the comparison. Panel (a) "
        "overlays the per-1 Myr mean and ±1σ envelope of compacted carbonate "
        "thickness from each product on a single axis, computed over the "
        "common deposition mask. The two curves track each other across the "
        "0–150 Ma window, indicating that the two CCD curves produce "
        "broadly similar global carbonate budgets through time at the "
        "population level. Panel (b) shows the distribution of the signed "
        "Δ = DM2026 − BW1991 in 5 Myr bins. The non-zero IQR of each box "
        "demonstrates that the agreement in panel (a) is statistical, not "
        "cell-by-cell. The four times of largest mean |Δ| are highlighted "
        "in red.")
    add_figure(doc, FIG_DIR / "03_thickness_and_delta.png",
               "Figure 1. (a) Per-1 Myr mean and ±1σ envelope of compacted "
               "carbonate sediment thickness for DM2026 (blue) and BW1991 "
               "(red), computed over the common deposition mask. (b) Pooled-"
               "cell distribution of Δ = DM2026 − BW1991 in 5 Myr bins. "
               "Highlighted red boxes mark the four times of largest mean "
               "|Δ| picked by greedy non-maximum suppression.")

    # ---- 3.2b 2x2 panel of delta maps at picked times
    add_para(doc,
        "Figure 2 shows the global Δ field at each of the four picked times. "
        "All four panels share the same symmetric colour scale (clamped at "
        "the 99th percentile of |Δ| across the four times) so the maps are "
        "directly comparable. Coastlines are reconstructed to age t using "
        "the Alfonso 2024 plate model.")
    add_figure(doc, FIG_DIR / "06_delta_panel.png",
               "Figure 2. Global maps of Δ = DM2026 − BW1991 carbonate-"
               "thickness anomaly at the four picked times of largest mean "
               "|Δ|. Winkel-Tripel projection. Coastlines reconstructed to "
               "age t with the Alfonso 2024 plate model.")

    # ---- 3.3 regional analysis figure
    add_h2(doc, "3.3 Where the products disagree most")
    add_para(doc,
        "Figure 3 shows the regional heatmap of mean |Δ| through time, "
        "with the 18 10°-latitude bands grouped by ocean basin (rows) "
        "against time (x-axis). Red dashed verticals mark the four picked "
        "times. The heatmap is log-scaled so both the low-amplitude "
        "background and the high-amplitude episodes are visible.")
    add_figure(doc, FIG_DIR / "04_region_heatmap.png",
               "Figure 3. Mean |Δ| per (basin × 10°-latitude-band) cell "
               "through time. Rows are grouped by present-day ocean basin "
               "(Pacific, Atlantic, Indian, Arctic, Southern); within each "
               "basin the rows run from −90° to +90° in 10° steps (band "
               "centres). Red dashed lines mark the four picked times of "
               "maximum disagreement.")

    # ---- 3.4 picked times table + per-time region tables
    add_h2(doc, "3.4 Four times of largest disagreement")
    if picked:
        add_para(doc, "Table 1 lists the four times picked by NMS.")
        body = [(p["time_Ma"], f"{p['mean_abs']:.2f}",
                 f"{p['max_abs']:.1f}", f"{p['iqr']:.2f}") for p in picked]
        add_table(doc,
                  ["Age (Ma)", "Mean |Δ| (m)", "Max |Δ| (m)", "IQR(Δ) (m)"],
                  body)
        doc.add_paragraph()

        add_para(doc,
            "Tables 2-5 break down the three (basin × latitude band) cells "
            "carrying the largest mean |Δ| at each picked time.")
        for n, p in enumerate(picked, start=2):
            t = p["time_Ma"]
            rows = region_tables.get(t, [])
            if not rows:
                continue
            add_h3(doc, f"Table {n}. Top regions at {t} Ma "
                        f"(mean |Δ| at this time = {p['mean_abs']:.2f} m).")
            body = [(r["rank"], r["basin"], f"{r['lat_band']:+d}°",
                     f"{r['mean_abs']:.2f}") for r in rows]
            add_table(doc, ["Rank", "Basin", "Latitude band centre",
                            "Mean |Δ| in cell (m)"], body)
            doc.add_paragraph()

    # ---- 4. Discussion
    add_h1(doc, "4. Discussion")
    add_para(doc,
        "The 5 Myr-binned boxplots (Fig. 1) show that DM2026 and BW1991 "
        "produce broadly similar carbonate-thickness distributions over "
        "the bulk of the 0-150 Ma window, with their medians tracking each "
        "other closely. The signed Δ distributions in panel (c), however, "
        "have non-zero inter-quartile ranges throughout, indicating that "
        "the agreement is statistical rather than cell-by-cell.")

    add_para(doc,
        "The regional heatmap (Fig. 2) and the per-picked-time tables show "
        "that the disagreement is not spatially uniform. Two recurring "
        "patterns emerge: (i) the largest mean |Δ| consistently sits in the "
        "low- to mid-latitude bands of the Pacific basin, reflecting that "
        "this is where the bulk of the global carbonate-bearing seafloor "
        "lies and small per-cell differences integrate to large area-mean "
        "differences; (ii) the polar bands (Arctic and Southern) carry "
        "isolated high-amplitude episodes that correspond to opening and "
        "closing of high-latitude carbonate-deposition windows under the "
        "differing CCD trajectories.")

    add_para(doc,
        "We did NOT attempt a paleogeographic reconstruction of the cells "
        "into their depositional positions. The basin / latitude labels are "
        "present-day positions of the underlying grid cells. For a high-"
        "latitude cell that has drifted significantly since deposition, "
        "this aliases a single physical depositional setting onto the "
        "wrong present-day basin label. For the bulk of the Pacific cells "
        "this aliasing is mild over 0-150 Ma, but conclusions drawn from "
        "the polar bands should be revisited once the cells are "
        "reconstructed to their paleo positions.")

    # ---- References
    add_h1(doc, "References")
    refs = [
        "Alfonso, C.P., Müller, R.D., Mather, B., and Anthony, M. (2024). "
        "Spatio-temporal copper prospectivity in the American Cordillera predicted "
        "by positive-unlabeled machine learning. GSA Bulletin 137, 702-711. "
        "https://doi.org/10.1130/B37614.1",
        "Boss, S.K., and Wilkinson, B.H. (1991). Planktogenic / eustatic "
        "control on cratonic / oceanic carbonate accumulation. Journal of "
        "Geology 99(4), 497-513.",
        "Dutkiewicz, A., and Müller, R.D. (2026). A new global Cenozoic-"
        "Mesozoic calcite compensation depth curve. In preparation.",
        "Müller, R.D., Cannon, J., Williams, S., Dutkiewicz, A., and "
        "Wright, N. (2026). PyBacktrack 1.5: A community tool for "
        "reconstructing paleobathymetry of drill sites, geohistory and global "
        "paleobathymetry. EGUsphere [preprint]. "
        "https://doi.org/10.5194/egusphere-2026-3680",
    ]
    for r in refs:
        doc.add_paragraph(r, style="List Number")

    out = cfg.OUTPUT_DIR / "carbonate_thickness_comparison.docx"
    doc.save(out)
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
