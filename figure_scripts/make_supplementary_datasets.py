#!/usr/bin/env python3
"""
Build the two supplementary Excel datasets from the current workflow outputs.

    Supplementary_Dataset_S1_CCD.xlsx           global + regional CCD, basin areas
    Supplementary_Dataset_S2_carbon_model.xlsx  deep-Earth carbon-cycle model results

Both were previously assembled by hand, which is how they drifted out of step with
the analysis. Regenerate them with this script after any rerun:

    python make_supplementary_datasets.py

Every sheet records where its numbers came from, and each workbook carries a
Provenance sheet giving the source file and its modification time. Sources that
depend on the CCD are checked against the current CCD curve and reported as STALE
if they predate it, so a dataset built on a half-finished rerun says so instead of
looking fine. Sources that do not depend on the CCD (the reconstructed basin areas
and area fractions, which are plate-model geometry) are labelled static and are not
age-checked; flagging those too would only teach the reader to ignore the column.
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

# ---- path resolution (same convention as the other Paper/ scripts) ----------
_HERE = Path(__file__).resolve().parent


def _workflow_root() -> Path:
    for _p in (_HERE.parent / "CCD_workflow_clean", _HERE.parent, _HERE):
        if (_p / "steps").is_dir() and (_p / "ccdworkflow").is_dir():
            return _p
    raise SystemExit("cannot locate the CCD workflow root (expected steps/ + ccdworkflow/)")


CW = _workflow_root()
FIGDIR = _HERE / "Figures" if (_HERE / "Figures").is_dir() else CW / "figures"
OUT = _HERE if (_HERE / "Figures").is_dir() else CW

HYBRID = FIGDIR / "CCD_hybrid_DM2026.txt"
STEP1 = CW / "steps/step1_timescale_conversion/outputs"
XLSX = CW / "data/ccds_and_oceanbasin_fractions.xlsx"
AREAS = (CW / "steps/step2_ocean_basin_areas/ocean_basin_area_outputs"
            / "ocean_basin_areas_0_70Ma_million_km2.csv")
NB = CW / "steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Outputs"
ATM = NB / "Notebook05/csv/05_atmospheric_influx_all_sources.csv"
PLATE = NB / "Notebook02/csv/02_plate_influx.csv"
NET = NB / "Notebook05/csv/05_net_carbon_outflux.csv"

OBS_MAX_MA = 52.0            # end of the observational segment of the global CCD

# ---- shared look, matching the published workbooks -------------------------
BODY = Font(name="Arial", size=10)
HEAD = Font(name="Arial", size=10, bold=True)
TITLE = Font(name="Arial", size=13, bold=True)
FILL = PatternFill("solid", fgColor="DDE6F0")
CENTRE = Alignment(horizontal="center", vertical="center")

# (sheet, source, depends_on_ccd)
_provenance: list[tuple[str, Path, bool]] = []


def _record(sheet: str, src: Path, ccd_dependent: bool = True) -> Path:
    _provenance.append((sheet, src, ccd_dependent))
    return src


def _need(*paths: Path) -> None:
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise SystemExit("missing input(s):\n  " + "\n  ".join(str(p) for p in missing))


def _mtime(p: Path) -> _dt.datetime:
    return _dt.datetime.fromtimestamp(p.stat().st_mtime)


def read_xy(path: Path, names) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", engine="python", comment="#",
                       header=None, names=names).dropna()


# ---- sheet writers ---------------------------------------------------------
def write_readme(ws, title: str, lines) -> None:
    ws["A1"] = title
    ws["A1"].font = TITLE
    for i, line in enumerate(lines, start=3):
        ws.cell(row=i, column=1, value=line).font = BODY
    ws.column_dimensions["A"].width = 118


def write_table(ws, header, rows, widths, formats, groups=None) -> None:
    """One data sheet. `groups` gives (label, span) pairs for a two-row merged
    header, as used by the carbon-model sheets."""
    if groups:
        ws.cell(row=1, column=1, value=header[0])
        ws.merge_cells(start_row=1, start_column=1, end_row=2, end_column=1)
        col = 2
        for label, span in groups:
            ws.cell(row=1, column=col, value=label)
            ws.merge_cells(start_row=1, start_column=col, end_row=1, end_column=col + span - 1)
            for k in range(span):
                ws.cell(row=2, column=col + k, value=("min", "mean", "max")[k])
            col += span
        head_rows, first_data = (1, 2), 3
        ws.freeze_panes = "B3"
    else:
        for j, h in enumerate(header, start=1):
            ws.cell(row=1, column=j, value=h)
        head_rows, first_data = (1,), 2
        ws.freeze_panes = "A2"

    for r in head_rows:
        for j in range(1, len(widths) + 1):
            c = ws.cell(row=r, column=j)
            c.font, c.fill, c.alignment = HEAD, FILL, CENTRE

    for i, row in enumerate(rows, start=first_data):
        for j, v in enumerate(row, start=1):
            c = ws.cell(row=i, column=j, value=v)
            c.font = BODY
            c.number_format = formats[j - 1]

    for j, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(j)].width = w


def write_provenance(ws, reference: Path) -> None:
    ws["A1"] = "Provenance of this workbook"
    ws["A1"].font = TITLE
    ws["A3"] = (f"Generated {_dt.datetime.now():%Y-%m-%d %H:%M} by Paper/"
                "make_supplementary_datasets.py. Rerun that script after any change to "
                "the workflow so these numbers stay in step with the analysis.")
    ws["A3"].font = BODY
    header = ["Sheet", "Source file", "Source last modified", "Status"]
    ref_t = _mtime(reference)
    rows = []
    for sheet, src, ccd_dependent in _provenance:
        t = _mtime(src)
        if not ccd_dependent:
            status = "static input, independent of the CCD"
        elif t >= ref_t - _dt.timedelta(minutes=5):
            status = "current"
        else:
            status = f"STALE - predates {reference.name}, rerun the step that writes it"
        rows.append([sheet, str(src.relative_to(CW.parent)), t.strftime("%Y-%m-%d %H:%M"), status])
    for j, h in enumerate(header, start=1):
        c = ws.cell(row=5, column=j, value=h)
        c.font, c.fill, c.alignment = HEAD, FILL, CENTRE
    for i, row in enumerate(rows, start=6):
        for j, v in enumerate(row, start=1):
            ws.cell(row=i, column=j, value=v).font = BODY
    for j, w in enumerate((28, 74, 22, 52), start=1):
        ws.column_dimensions[get_column_letter(j)].width = w
    return [r for r in rows if r[3].startswith("STALE")]


# ---- Dataset S1 ------------------------------------------------------------
def build_s1() -> Path:
    _need(HYBRID, XLSX, AREAS,
          STEP1 / "ATL_CCD_GTS2020_1my.txt", STEP1 / "PAC_CCD_GTS2020_1my.txt",
          STEP1 / "IND_CCD_GTS2020_1my.txt")
    wb = Workbook()

    write_readme(wb.active, "Supplementary Dataset S1 — Global and regional carbonate "
                            "compensation depth (CCD) and ocean-basin areas", [
        "Dutkiewicz & Müller — Cenozoic global CCD (in prep).",
        "",
        "Contents:",
        "  • Global CCD — area-weighted observational global CCD (0–52 Ma) extended to 170 Ma as a sea-level-",
        "    calibrated, planktogenic-corrected hybrid estimate, with lower/upper uncertainty bounds.",
        "  • Regional CCDs — the Atlantic, Pacific and Indian Ocean CCD reconstructions combined into the global curve.",
        "  • Ocean basin area fractions — time-dependent fractional abyssal area used to area-weight the regional curves.",
        "    A basin's fraction is blank once its CCD record ends: the fractions are renormalised over the basins",
        "    that have data at that age, so the Indian Ocean area is split evenly between the two sampled basins",
        "    before 24 Ma, which is algebraically identical to imputing the unsampled Indian CCD as the midpoint",
        "    of the Atlantic and Pacific curves.",
        "  • Ocean basin areas — absolute reconstructed abyssal areas (10^6 km^2).",
        "",
        "Units and conventions:",
        "  • CCD in metres; negative = depth below present sea level (a deeper CCD is a more negative number).",
        "  • Area fractions are dimensionless (sum to 1 across the three basins at each age).",
        "  • Ages in Ma (millions of years before present).",
        "  • All CCD curves are on the GTS2020 timescale (Ogg, 2020). The published Atlantic curve is on GTS2012",
        "    and the 36–52 Ma Pacific segment on Cande & Kent (1995); both are converted by step 1 of the workflow.",
        "",
        "Sources (files in the CCD_workflow_clean repository):",
        "  • Global CCD: Paper/Figures/CCD_hybrid_DM2026.txt (0–52 Ma observational; 52–170 Ma hybrid).",
        "  • Regional CCDs: steps/step1_timescale_conversion/outputs/{ATL,PAC,IND}_CCD_GTS2020_1my.txt",
        "  • Area fractions: data/ccds_and_oceanbasin_fractions.xlsx",
        "  • Ocean areas: steps/step2_ocean_basin_areas/ocean_basin_area_outputs/",
        "    ocean_basin_areas_0_70Ma_million_km2.csv",
        "",
        "See the Provenance sheet for the exact files used and when they were last written.",
    ])
    wb.active.title = "Read me"

    # Global CCD
    h = read_xy(_record("Global CCD", HYBRID), ["age", "ccd", "lo", "hi"])
    rows = [[float(a), float(c), float(lo), float(hi),
             "observational (area-weighted)" if a <= OBS_MAX_MA else "hybrid (sea-level-calibrated)"]
            for a, c, lo, hi in h.itertuples(index=False)]
    write_table(wb.create_sheet("Global CCD"),
                ["Age (Ma)", "Global CCD (m)", "Lower bound (m)", "Upper bound (m)", "Segment"],
                rows, [10, 16, 16, 16, 32], ["0.000", "0.0", "0.0", "0.0", "General"])

    # Regional CCDs, on the common 1 Myr grid written by step 1
    reg = {}
    for key, fn in (("Atlantic", "ATL_CCD_GTS2020_1my.txt"),
                    ("Pacific", "PAC_CCD_GTS2020_1my.txt"),
                    ("Indian", "IND_CCD_GTS2020_1my.txt")):
        df = read_xy(_record(f"Regional CCDs ({key})", STEP1 / fn), ["age", "ccd"])
        reg[key] = dict(zip(df["age"].round(3), df["ccd"]))
    ages = sorted(set().union(*(set(v) for v in reg.values())))
    rows = [[a] + [reg[k].get(a) for k in ("Atlantic", "Pacific", "Indian")] for a in ages]
    write_table(wb.create_sheet("Regional CCDs"),
                ["Age (Ma)", "Atlantic CCD (m)", "Pacific CCD (m)", "Indian CCD (m)"],
                rows, [10, 18, 18, 18], ["0.000", "0.0", "0.0", "0.0"])

    # Area fractions (columns 6-9 of the spreadsheet: master age grid, then ATL/PAC/IND)
    df = pd.read_excel(_record("Ocean basin area fractions", XLSX, ccd_dependent=False), header=1)
    f = df.iloc[:, 6:10].apply(pd.to_numeric, errors="coerce")
    f = f[f.iloc[:, 0].notna()]
    # A basin's cell is blank once its CCD record ends, because the fractions are
    # renormalised over the basins that have data at that age. Keep the blanks: they
    # are the weighting, not missing values, and dropping those rows would discard
    # every age older than the Indian Ocean record.
    rows = [[float(a)] + [None if pd.isna(v) else float(v) for v in (x, y, z)]
            for a, x, y, z in f.itertuples(index=False)]
    write_table(wb.create_sheet("Ocean basin area fractions"),
                ["Age (Ma)", "Atlantic fraction", "Pacific fraction", "Indian fraction"],
                rows, [10, 18, 18, 18], ["0.000"] * 4)

    # Absolute areas
    a = pd.read_csv(_record("Ocean basin areas", AREAS, ccd_dependent=False))
    a.columns = [c.strip() for c in a.columns]
    col = {c.split(" (")[0]: c for c in a.columns}
    rows = [[float(r[col["Time"]]), float(r[col["Atlantic"]]), float(r[col["Pacific"]]),
             float(r[col["Indian"]]), float(r[col["Sum"]])] for _, r in a.iterrows()]
    write_table(wb.create_sheet("Ocean basin areas"),
                ["Time (Ma)", "Atlantic (10^6 km2)", "Pacific (10^6 km2)",
                 "Indian (10^6 km2)", "Sum (10^6 km2)"],
                rows, [10, 20, 20, 20, 20], ["0.0"] * 5)

    stale = write_provenance(wb.create_sheet("Provenance"), HYBRID)
    path = OUT / "Supplementary_Dataset_S1_CCD.xlsx"
    wb.save(path)
    return path, stale


# ---- Dataset S2 ------------------------------------------------------------
def _triplet(df, stem):
    return [df[f"{stem}_{q}"].to_numpy(float) for q in ("min", "mean", "max")]


def _mi_triplet(df, stem):
    return [df[(stem, q)].to_numpy(float) for q in ("min", "mean", "max")]


def build_s2() -> Path:
    _need(ATM, PLATE, NET)
    wb = Workbook()

    write_readme(wb.active, "Supplementary Dataset S2 — Deep-Earth carbon-cycle model results", [
        "Dutkiewicz & Müller — recomputed by driving the Müller et al. (2024) degassing model with the new global CCD",
        "on the Alfonso et al. (2024) plate model.",
        "",
        "Contents:",
        "  • Atmospheric outflux by source — solid-Earth CO2 degassing to the atmosphere, separated by reservoir,",
        "    plus the gross atmospheric outflux.",
        "  • Plate influx by reservoir — carbon carried into the plate / sequestered, separated by reservoir,",
        "    plus the gross (total) plate influx.",
        "  • Net atmospheric flux — gross atmospheric outflux minus gross plate influx, in two versions:",
        "    excluding and including the growing pelagic (deep-sea) carbonate-sediment reservoir in the storage term.",
        "",
        "Units and conventions:",
        "  • All fluxes in Mt C yr^-1 (megatonnes of carbon per year).",
        "  • Atmospheric outflux positive = CO2 released to the atmosphere (degassing).",
        "  • Plate influx positive = carbon carried into / stored in the plate.",
        "  • Net atmospheric flux positive = net CO2 to the atmosphere.",
        "  • min / mean / max span the model's parameter-uncertainty envelope.",
        "  • Rift outflux uses the 'biased-rift' formulation (consistent with the gross outflux reported here).",
        "  • Serpentinite = total (bending-related + mid-ocean-ridge) serpentinite carbon.",
        "",
        "Sources (files in the CCD_workflow_clean repository, steps/step10_carbon_cycle_degassing/",
        "Alfonso_etal_2024_DM26/Outputs):",
        "  • Atmospheric outflux: Notebook05/csv/05_atmospheric_influx_all_sources.csv",
        "  • Plate influx: Notebook02/csv/02_plate_influx.csv",
        "  • Net flux: Notebook05/csv/05_net_carbon_outflux.csv",
        "",
        "See the Provenance sheet for the exact files used and when they were last written.",
    ])
    wb.active.title = "Read me"

    atm = pd.read_csv(_record("Atmospheric outflux", ATM)).rename(columns={"Age (Ma)": "age"})
    ages = atm["age"].to_numpy(float)
    groups = [("Mid-ocean ridge", "ridge_outflux"), ("Arc (subduction)", "subduction_outflux"),
              ("Carbonate platform", "carbonate_platform_outflux"), ("Rift", "rift_outflux_biased"),
              ("Intraplate", "intraplate_volcanism_outflux"),
              ("Gross atmospheric outflux", "gross_atmospheric_outflux_biased_rift")]
    cols = [ages] + [c for _, stem in groups for c in _triplet(atm, stem)]
    write_table(wb.create_sheet("Atmospheric outflux"), ["Age (Ma)"],
                list(map(list, zip(*cols))), [10] + [11] * (len(cols) - 1),
                ["General"] + ["0.000"] * (len(cols) - 1),
                groups=[(label, 3) for label, _ in groups])

    pl = pd.read_csv(_record("Plate influx", PLATE), index_col=0, header=[0, 1])
    pl.index = pd.to_numeric(pl.index, errors="coerce")
    pl = pl[np.isfinite(pl.index)]
    pgroups = [("Pelagic carbonate sediment", "sediments"), ("Altered oceanic crust", "crust"),
               ("Serpentinite (total)", "serpentinite_total"),
               ("Organic sediments", "organic_sediments"), ("Gross plate influx", "total_influx")]
    cols = [pl.index.to_numpy(float)] + [c for _, stem in pgroups for c in _mi_triplet(pl, stem)]
    write_table(wb.create_sheet("Plate influx"), ["Age (Ma)"],
                list(map(list, zip(*cols))), [10] + [11] * (len(cols) - 1),
                ["General"] + ["0.000"] * (len(cols) - 1),
                groups=[(label, 3) for label, _ in pgroups])

    net = pd.read_csv(_record("Net atmospheric flux", NET), index_col=0, header=[0, 1])
    net.index = pd.to_numeric(net.index, errors="coerce")
    net = net[np.isfinite(net.index)]
    ngroups = [("Net (excl. pelagic sediment reservoir)", "net_carbon_outflux"),
               ("Net (incl. pelagic sediment reservoir)",
                "net_carbon_outflux_with_sediment_storage_biased")]
    cols = [net.index.to_numpy(float)] + [c for _, stem in ngroups for c in _mi_triplet(net, stem)]
    write_table(wb.create_sheet("Net atmospheric flux"), ["Age (Ma)"],
                list(map(list, zip(*cols))), [10] + [11] * (len(cols) - 1),
                ["General"] + ["0.000"] * (len(cols) - 1),
                groups=[(label, 3) for label, _ in ngroups])

    stale = write_provenance(wb.create_sheet("Provenance"), HYBRID)
    path = OUT / "Supplementary_Dataset_S2_carbon_model.xlsx"
    wb.save(path)
    return path, stale


def main() -> None:
    global _provenance
    for build in (build_s1, build_s2):
        _provenance = []
        path, stale = build()
        print(f"wrote {path}")
        for sheet, src, when, status in stale:
            print(f"  !! {sheet}: {status}  ({src}, {when})")
    print("\nThe reference for 'current' is Paper/Figures/CCD_hybrid_DM2026.txt.")


if __name__ == "__main__":
    main()
