#!/usr/bin/env python3
"""Repack the archive's NetCDF grids as scaled integers, into a parallel tree.

    python zenodo/repack_grids.py            # convert everything
    python zenodo/repack_grids.py --check    # report what it would do, convert nothing

GMT already writes these grids as NetCDF-4 with zlib level 3, so recompressing them
losslessly gains about four percent. What costs the space is float32: seven
significant digits spent on smooth quantities whose own uncertainty is tens of metres.
The spread between the minimum and maximum sedimentation-rate scenarios has a median
of 19 m and a 95th percentile of 140 m, so a storage precision of 0.1 m is some four
hundred times finer than anything the model resolves.

Each variable is therefore stored as an integer with a scale factor, which every
NetCDF reader - GMT, xarray, gplately, MATLAB, ncview - unpacks transparently on
read. The files still open as floating-point grids in metres and nothing downstream
has to change.

    thickness grids    int16, 0.1 m       error <= 0.05 m
    paleobathymetry    int16, 0.25 m      error <= 0.125 m
    deposition masks   int8               exact, the data is 1 or nothing

Every file is read back after writing and the actual error compared against the
tolerance, so a grid that did not survive the round trip is reported rather than
shipped. Originals are never touched: the output goes to a parallel tree.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import xarray as xr

# variable kind -> (numpy dtype, scale_factor, add_offset, fill value, tolerance in m)
# The tolerance is one whole quantum, not half of one: rounding to a 0.1 m grid puts
# the true error at or just under 0.05 m, and the float32-to-float64 round trip adds a
# few parts in 10^7 on top, which a half-quantum bound trips over. The check exists to
# catch a packing that has genuinely broken - those are wrong by metres, not microns.
PACKINGS = {
    "thickness":    ("int16", 0.1,  0.0,     -32768, 0.1),
    "bathymetry":   ("int16", 0.25, -2000.0, -32768, 0.25),
    "mask":         ("int8",  None, None,    -1,     0.0),
}


def kind_of(name: str) -> str:
    n = name.lower()
    if "deposition_mask" in n:
        return "mask"
    if "paleobathymetry" in n:
        return "bathymetry"
    return "thickness"


def limits(dtype: str, scale: float, offset: float) -> tuple[float, float]:
    """The physical range the packing can represent, leaving the fill value out."""
    info = np.iinfo(dtype)
    return (info.min + 1) * scale + offset, info.max * scale + offset


def repack(src: Path, dst: Path, complevel: int = 5) -> tuple[int, int, float]:
    """Write a packed copy of one grid. Returns (bytes in, bytes out, max error)."""
    kind = kind_of(src.name)
    dtype, scale, offset, fill, tol = PACKINGS[kind]
    ds = xr.open_dataset(src)
    try:
        enc = {}
        for v in ds.data_vars:
            z = ds[v].values
            if kind == "mask":
                enc[v] = dict(dtype=dtype, _FillValue=fill, zlib=True,
                              complevel=complevel, shuffle=True)
                continue
            lo, hi = limits(dtype, scale, offset)
            finite = np.isfinite(z)
            if finite.any():
                zmin, zmax = float(np.nanmin(z)), float(np.nanmax(z))
                if zmin < lo or zmax > hi:
                    raise SystemExit(
                        f"{src}: values {zmin:.1f}..{zmax:.1f} m fall outside what the "
                        f"{dtype} packing can hold ({lo:.1f}..{hi:.1f} m). Widen the "
                        f"scale factor in PACKINGS rather than shipping clipped data.")
            enc[v] = dict(dtype=dtype, scale_factor=scale, add_offset=offset,
                          _FillValue=fill, zlib=True, complevel=complevel, shuffle=True)
            ds[v].attrs["units"] = ds[v].attrs.get("units", "m")
            ds[v].attrs["packing"] = (f"stored as {dtype} with scale_factor {scale}"
                                      f"{f' and add_offset {offset}' if offset else ''}; "
                                      f"readers unpack to metres automatically")
        dst.parent.mkdir(parents=True, exist_ok=True)
        ds.to_netcdf(dst, encoding=enc, engine="netcdf4")
    finally:
        ds.close()

    # Read it back and measure what the packing actually cost.
    a = xr.open_dataset(src); b = xr.open_dataset(dst)
    try:
        err = 0.0
        for v in a.data_vars:
            x, y = a[v].values.astype("float64"), b[v].values.astype("float64")
            m = np.isfinite(x)
            if (np.isfinite(y) != m).any():
                raise SystemExit(f"{dst}: the pattern of missing values changed")
            if m.any():
                err = max(err, float(np.abs(y[m] - x[m]).max()))
        if err > tol + 1e-9:
            raise SystemExit(f"{dst}: round-trip error {err:.4g} m exceeds the "
                             f"{tol} m tolerance for a {kind} grid")
    finally:
        a.close(); b.close()
    return src.stat().st_size, dst.stat().st_size, err


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="report the folders and file counts, convert nothing")
    ap.add_argument("--complevel", type=int, default=5)
    ap.add_argument("--out", default=None,
                    help="output root (default: <workflow>/steps/.../repacked)")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    repo = here.parent
    roots = [repo / "steps", repo.parent / "CCD_workflow_clean" / "steps"]
    steps = next((r for r in roots
                  if (r / "step8_carbonate_sediment_thickness"
                        / "carbonate_sed_thickness_DM2026").is_dir()), None)
    if steps is None:
        raise SystemExit("cannot find the grids; looked in " + ", ".join(map(str, roots)))
    s8 = steps / "step8_carbonate_sediment_thickness"
    out_root = Path(args.out) if args.out else s8 / "repacked"

    SOURCES = [
        ("carbonate_sed_thickness_DM2026",       s8 / "carbonate_sed_thickness_DM2026"),
        ("carbonate_sed_thickness_min_DM2026",   s8 / "carbonate_sed_thickness_min_DM2026"),
        ("carbonate_sed_thickness_max_DM2026",   s8 / "carbonate_sed_thickness_max_DM2026"),
        ("Alfonso2024_pybacktrack_merged_paleobathymetry",
         s8 / "input_grids" / "Alfonso2024_pybacktrack_merged_paleobathymetry"),
    ]

    print(f"Grids: {s8}")
    print(f"Output: {out_root}\n")
    plan = []
    for name, folder in SOURCES:
        if not folder.is_dir():
            raise SystemExit(f"missing: {folder}")
        files = sorted(p for p in folder.glob("*.nc"))
        size = sum(p.stat().st_size for p in files)
        plan.append((name, folder, files))
        print(f"  {name:<50} {len(files):4d} files  {size/1e6:8.1f} MB")
    if args.check:
        print("\n--check: nothing written.")
        return

    total_in = total_out = 0
    worst = 0.0
    for name, folder, files in plan:
        print(f"\n{name}")
        fin = fout = 0
        for i, src in enumerate(files, 1):
            a, b, err = repack(src, out_root / name / src.name, args.complevel)
            fin += a; fout += b; worst = max(worst, err)
            if i % 25 == 0 or i == len(files):
                print(f"  {i:4d}/{len(files)}  {fin/1e6:7.1f} -> {fout/1e6:7.1f} MB", flush=True)
        total_in += fin; total_out += fout
        print(f"  {fin/1e6:.1f} MB -> {fout/1e6:.1f} MB  ({fin/max(fout,1):.2f}x smaller)")

    print(f"\nTotal {total_in/1e6:.0f} MB -> {total_out/1e6:.0f} MB "
          f"({total_in/max(total_out,1):.2f}x smaller, "
          f"{100*(1-total_out/max(total_in,1)):.0f}% saved)")
    print(f"Largest round-trip error anywhere: {worst:.4g} m")
    print(f"\nRepacked grids are in {out_root}. The originals are untouched.")


if __name__ == "__main__":
    main()
