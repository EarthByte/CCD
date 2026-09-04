#!/usr/bin/env python3
"""
Step 1 — normalise every input curve to the GTS2020 timescale.

The regional CCD curves and the sea-level records were published on three
different geomagnetic polarity timescales. Combining them without conversion
mixes age models, so every series is first mapped onto GTS2020 (Ogg, 2020) and
then resampled onto a common grid.

Timescale of each input
-----------------------
    Atlantic CCD      Dutkiewicz & Muller (2022)      GTS2012      -> convert
    Pacific CCD       Dalvand et al. (2025), 0-35 Ma  GTS2020      -> as is
                      Palike et al. (2012), 36-52 Ma  Cande & Kent 1995 -> convert
    Indian CCD        Dalvand et al. (2025)           GTS2020      -> as is
    Sea level         Miller et al. (2024), 0-66.6 Ma GTS2020      -> as is
                      Haq et al., older than 66.6 Ma  GTS2012      -> convert
    Haq short-term    Haq et al. (1987)               GTS2012      -> convert
    Haq long-term     Haq et al. (1987)               GTS2012      -> convert

The two Haq records are the validation pair used in step 4 (the envelope operator
is applied to the short-term curve and checked against the published long-term
curve), so both are converted, keeping them in a common frame.

Only the segments that need it are converted; a segment already on GTS2020 is
passed through untouched, so no spurious interpolation is introduced.

The age mapping is done with the EarthByte `timescale_conversion.py` tool and its
`timescales.txt` chron table, both copied into this folder so the workflow is
self-contained (originals: scripts/1_timescale_conversion/timescale_conversion/).
Ages are mapped by piecewise-linear interpolation between the chron boundaries
common to the two timescales.

Resampling
----------
An age remap leaves the samples unevenly spaced, so each series is put back onto
a uniform grid afterwards: 1 Myr for the CCD curves, 0.1 Myr for the sea-level
series (the rolling-quantile envelope in step 4 needs the dense sampling and
would be destroyed by a 1 Myr grid).

Inputs:
    data/ccds_and_oceanbasin_fractions.xlsx      regional CCD columns
    data/sealevel/Miller_Haq_SeaLevel_ShortTerm_hybrid.tsv
    data/sealevel/Haq_SeaLevel_ShortTerm_hybrid.tsv
    data/sealevel/Haq87_Longterm_v3.txt

Outputs (outputs/):
    ATL_CCD_GTS2020_1my.txt        Age_Ma, CCD_m
    PAC_CCD_GTS2020_1my.txt
    IND_CCD_GTS2020_1my.txt
    sealevel_shortterm_hybrid_GTS2020.txt   Age_Ma, Sea_level_m
    Haq_shortterm_hybrid_GTS2020.txt
    Haq87_longterm_GTS2020.txt
    timescale_conversion_report.txt         what moved, and by how much
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from ccdworkflow import config
from ccdworkflow.io import write_table

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
TIMESCALES_TXT = HERE / "timescales.txt"

# names as spelled in timescale_conversion.py
GTS2020 = "Ogg_2020"
GTS2012 = "Ogg_2012"
CK95 = "Cande_Kent_1995"

# the Miller / Haq splice in the hybrid short-term sea-level record
SL_SPLICE_MA = 66.6
# the Dalvand / Palike splice in the stitched equatorial Pacific CCD
PAC_SPLICE_MA = 35.5


def _load_tool():
    """Import the EarthByte timescale_conversion.py sitting next to this script."""
    spec = importlib.util.spec_from_file_location(
        "earthbyte_timescale_conversion", HERE / "timescale_conversion.py"
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _read_timescales(tool):
    """read_timescales() drops a timescales.log in the working directory, so run it
    from outputs/ rather than wherever the workflow happens to be invoked."""
    OUT.mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    try:
        os.chdir(OUT)
        scales = tool.read_timescales(str(TIMESCALES_TXT))
    finally:
        os.chdir(cwd)
    return {name: scales[i] for i, name in enumerate(tool.timescale_names)}


def convert_ages(tool, scales, ages, source):
    """Map ages from `source` onto GTS2020. Ages beyond the chron table are left
    alone by the tool; they are counted and reported rather than silently dropped."""
    if source == GTS2020:
        return np.asarray(ages, float), 0
    old, new = scales[source], scales[GTS2020]
    out, moved = [], 0
    for a in np.asarray(ages, float):
        changed, na = tool.get_new_age(float(a), old, new)
        out.append(na)
        moved += int(changed)
    return np.asarray(out, float), moved


def convert_series(tool, scales, ages, values, segments, name, log):
    """Convert a curve whose segments sit on different timescales.

    `segments` is a list of (age_lo, age_hi, timescale) covering the series in the
    ORIGINAL age frame, youngest first; each is converted independently and the
    pieces concatenated.

    The splice age between two segments is defined in the ORIGINAL frame, so once the
    older segment is converted its youngest samples can cross back over the boundary
    and interleave with the younger record (the GTS2012 -> GTS2020 map moves late
    Cretaceous ages about 1.3 Myr younger, for instance). Interleaving two independent
    records would put a sawtooth into the spliced curve, so the splice is re-imposed in
    the CONVERTED frame: any sample of an older segment that lands at or below the
    highest converted age already accepted is dropped. The 0.1 / 1 Myr regrid then
    bridges the resulting short gap linearly.
    """
    ages = np.asarray(ages, float)
    values = np.asarray(values, float)
    order = np.argsort(ages)
    ages, values = ages[order], values[order]

    new_ages, new_vals = [], []
    ceiling = -np.inf          # highest converted age accepted so far
    for lo, hi, ts in segments:
        m = (ages >= lo) & (ages <= hi)
        if not m.any():
            continue
        ca, moved = convert_ages(tool, scales, ages[m], ts)
        cv = values[m]
        shift = ca - ages[m]
        log.append(
            f"  {name:<12s} {lo:6.1f}-{hi:6.1f} Ma  from {ts:<16s} "
            f"n={m.sum():4d}  ages moved: {moved:4d}  "
            f"shift mean {shift.mean():+7.3f} Myr, max |{np.abs(shift).max():.3f}| Myr"
        )
        if np.isfinite(ceiling):
            crossed = ca <= ceiling + 1e-9
            if crossed.any():
                log.append(
                    f"  {name:<12s} splice re-imposed at {ceiling:.3f} Ma (GTS2020): "
                    f"{int(crossed.sum())} converted sample(s) crossed back below it "
                    f"and were dropped; gap to the next sample "
                    f"{ca[~crossed].min() - ceiling:.3f} Myr"
                )
            ca, cv = ca[~crossed], cv[~crossed]
            if ca.size == 0:
                continue
        new_ages.append(ca)
        new_vals.append(cv)
        ceiling = max(ceiling, float(ca.max()))
    a = np.concatenate(new_ages)
    v = np.concatenate(new_vals)
    o = np.argsort(a)
    a, v = a[o], v[o]
    # a coarse chron mapping can also fold two neighbouring samples onto one age
    keep = np.concatenate([[True], np.diff(a) > 1e-9])
    dropped = int((~keep).sum())
    if dropped:
        log.append(f"  {name:<12s} dropped {dropped} sample(s) folded onto a duplicate age")
    return a[keep], v[keep]


def regrid(ages, values, step):
    """Put a converted series back onto a uniform grid, without extrapolating."""
    lo = np.ceil(ages.min() / step) * step
    hi = np.floor(ages.max() / step) * step
    grid = np.round(np.arange(lo, hi + step / 2, step), 6)
    return grid, np.interp(grid, ages, values)


def read_xlsx_regional():
    """Regional CCD curves as the synthesis uses them, from the shared spreadsheet
    (col 0/1 Atlantic, 2/3 Pacific, 4/5 Indian)."""
    df = pd.read_excel(config.CCD_FRACTIONS_XLSX, header=1)
    out = {}
    for name, ca, cv in (("Atlantic", 0, 1), ("Pacific", 2, 3), ("Indian", 4, 5)):
        a = pd.to_numeric(df.iloc[:, ca], errors="coerce")
        v = pd.to_numeric(df.iloc[:, cv], errors="coerce")
        m = np.isfinite(a) & np.isfinite(v)
        out[name] = (a[m].to_numpy(float), v[m].to_numpy(float))
    return out


def read_xy(path):
    a, v = [], []
    for line in open(path):
        if line.startswith(("#", ">")):
            continue
        s = line.split()
        if len(s) < 2:
            continue
        try:
            a.append(float(s[0]))
            v.append(float(s[1]))
        except ValueError:
            continue
    return np.asarray(a, float), np.asarray(v, float)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tool = _load_tool()
    scales = _read_timescales(tool)
    log = ["Timescale conversion to GTS2020 (Ogg, 2020)", ""]

    # ---- regional CCD curves -------------------------------------------------
    reg = read_xlsx_regional()
    plans = {
        "Atlantic": ([(0.0, 1e9, GTS2012)], "ATL_CCD_GTS2020_1my.txt"),
        "Pacific": ([(0.0, PAC_SPLICE_MA, GTS2020), (PAC_SPLICE_MA, 1e9, CK95)],
                    "PAC_CCD_GTS2020_1my.txt"),
        "Indian": ([(0.0, 1e9, GTS2020)], "IND_CCD_GTS2020_1my.txt"),
    }
    log.append("Regional CCD curves (resampled to 1 Myr after conversion):")
    for name, (segments, fname) in plans.items():
        a, v = reg[name]
        ca, cv = convert_series(tool, scales, a, v, segments, name, log)
        ga, gv = regrid(ca, cv, 1.0)
        write_table(OUT / fname,
                    pd.DataFrame({"Age_Ma": ga, "CCD_m": np.round(gv, 2)}),
                    header=["Age_Ma", "CCD_m"], float_format="%.2f")
        log.append(f"  {name:<12s} -> {fname}  ({len(ga)} rows, {ga.min():.0f}-{ga.max():.0f} Ma)")
    log.append("")

    # ---- sea level -----------------------------------------------------------
    log.append("Sea level (resampled to 0.1 Myr after conversion):")
    sl_in = config.SEALEVEL / "Miller_Haq_SeaLevel_ShortTerm_hybrid.tsv"
    a, v = read_xy(sl_in)
    ca, cv = convert_series(
        tool, scales, a, v,
        [(0.0, SL_SPLICE_MA, GTS2020), (SL_SPLICE_MA, 1e9, GTS2012)],
        "SL hybrid", log)
    ga, gv = regrid(ca, cv, 0.1)
    write_table(OUT / "sealevel_shortterm_hybrid_GTS2020.txt",
                pd.DataFrame({"Age_Ma": ga, "Sea_level_m": np.round(gv, 6)}),
                header=["Age_Ma", "Sea_level_m"], float_format="%.6f")
    log.append(f"  {'SL hybrid':<12s} -> sealevel_shortterm_hybrid_GTS2020.txt "
               f"({len(ga)} rows, {ga.min():.1f}-{ga.max():.1f} Ma)")

    a, v = read_xy(config.SEALEVEL / "Haq_SeaLevel_ShortTerm_hybrid.tsv")
    ca, cv = convert_series(tool, scales, a, v, [(0.0, 1e9, GTS2012)], "Haq ST", log)
    ga, gv = regrid(ca, cv, 0.1)
    write_table(OUT / "Haq_shortterm_hybrid_GTS2020.txt",
                pd.DataFrame({"Age_Ma": ga, "Sea_level_m": np.round(gv, 6)}),
                header=["Age_Ma", "Sea_level_m"], float_format="%.6f")
    log.append(f"  {'Haq ST':<12s} -> Haq_shortterm_hybrid_GTS2020.txt "
               f"({len(ga)} rows, {ga.min():.1f}-{ga.max():.1f} Ma)")

    a, v = read_xy(config.SEALEVEL / "Haq87_Longterm_v3.txt")
    ca, cv = convert_series(tool, scales, a, v, [(0.0, 1e9, GTS2012)], "Haq87 LT", log)
    ga, gv = regrid(ca, cv, 0.1)
    write_table(OUT / "Haq87_longterm_GTS2020.txt",
                pd.DataFrame({"Age_Ma": ga, "Sea_level_m": np.round(gv, 6)}),
                header=["Age_Ma", "Sea_level_m"], float_format="%.6f")
    log.append(f"  {'Haq87 LT':<12s} -> Haq87_longterm_GTS2020.txt "
               f"({len(ga)} rows, {ga.min():.1f}-{ga.max():.1f} Ma)")

    (OUT / "timescale_conversion_report.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))
    print(f"\n[step1] wrote {OUT}")


if __name__ == "__main__":
    main()
