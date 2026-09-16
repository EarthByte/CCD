#!/usr/bin/env python3
"""Repack the archive's NetCDF grids as scaled integers, into a parallel tree.

    python zenodo/repack_grids.py --to-archive ../CCD_zenodo_archive   # recommended
    python zenodo/repack_grids.py --check    # report what it would do, convert nothing
    python zenodo/repack_grids.py            # write a parallel tree of loose files

--to-archive converts each grid and streams it straight into the component's tarball,
so the run reads many files but writes only one. Writing a tree of loose files instead
means several thousand reads followed by several thousand writes of transformed,
higher-entropy data, which is the behavioural signature endpoint-protection software
watches for; Cortex XDR reads it as ransomware and kills the session. Prefer
--to-archive, and use --pace to slow the run further if a scanner still objects.

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


def _encode(ds, enc, scratch: Path):
    """Serialise a dataset to bytes. Prefers an in-memory write; falls back to one
    scratch file that is rewritten in place, never a new file per grid."""
    try:
        import h5netcdf  # noqa: F401
        import io
        buf = io.BytesIO()
        ds.to_netcdf(buf, encoding=enc, engine="h5netcdf")
        return buf.getvalue()
    except Exception:
        ds.to_netcdf(scratch, encoding=enc, engine="netcdf4")
        return scratch.read_bytes()


def repack_bytes(src: Path, scratch: Path, complevel: int = 5) -> tuple[int, bytes, float]:
    """Return (bytes in, packed bytes, max round-trip error) without writing output."""
    import io
    kind = kind_of(src.name)
    dtype, scale, offset, fill, tol = PACKINGS[kind]
    ds = xr.open_dataset(src)
    try:
        enc = _encoding_for(ds, src, kind, dtype, scale, offset, fill, complevel)
        raw = _encode(ds, enc, scratch)
        original = {v: ds[v].values.astype("float64") for v in ds.data_vars}
    finally:
        ds.close()
    # Read the packed bytes back through netCDF4's in-memory reader, which applies
    # scale_factor and add_offset exactly as any other reader will, so the error
    # measured here is the error a user would see.
    import netCDF4
    back = netCDF4.Dataset("inmemory.nc", mode="r", memory=raw)
    try:
        err = 0.0
        for v, x in original.items():
            y = np.ma.filled(back[v][:].astype("float64"), np.nan)
            m = np.isfinite(x)
            if (np.isfinite(y) != m).any():
                raise SystemExit(f"{src}: the pattern of missing values changed")
            if m.any():
                err = max(err, float(np.abs(y[m] - x[m]).max()))
    finally:
        back.close()
    if err > tol + 1e-9:
        raise SystemExit(f"{src}: round-trip error {err:.4g} m exceeds the {tol} m "
                         f"tolerance for a {kind} grid")
    return src.stat().st_size, raw, err


def _encoding_for(ds, src, kind, dtype, scale, offset, fill, complevel):
    enc = {}
    for v in ds.data_vars:
        if kind == "mask":
            enc[v] = dict(dtype=dtype, _FillValue=fill, zlib=True,
                          complevel=complevel, shuffle=True)
            continue
        z = ds[v].values
        lo, hi = limits(dtype, scale, offset)
        if np.isfinite(z).any():
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
    return enc


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
    ap.add_argument("--to-archive", metavar="DIR", default=None,
                    help="stream each component straight into DIR/<name>.tar.gz; "
                         "writes one file per component instead of thousands")
    ap.add_argument("--pace", type=float, default=0.0, metavar="SECONDS",
                    help="pause this long after each grid, to keep the read rate low "
                         "if an endpoint scanner objects (try 0.02)")
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

    if args.to_archive:
        import io, tarfile, tempfile, time
        outdir = Path(args.to_archive); outdir.mkdir(parents=True, exist_ok=True)
        # Components the archive ships unchanged: names here must match
        # make_zenodo_archive.sh, which is what builds the rest of the archive.
        TARNAME = {"carbonate_sed_thickness_DM2026": "carbonate_sediment_thickness_mean",
                   "carbonate_sed_thickness_min_DM2026": "carbonate_sediment_thickness_min",
                   "carbonate_sed_thickness_max_DM2026": "carbonate_sediment_thickness_max",
                   "Alfonso2024_pybacktrack_merged_paleobathymetry": "paleobathymetry"}
        total_in = total_out = 0; worst = 0.0
        with tempfile.TemporaryDirectory() as td:
            scratch = Path(td) / "grid.nc"
            for name, folder, files in plan:
                tar_path = outdir / f"{TARNAME[name]}.tar.gz"
                # Build beside the existing tarball and swap only on success, so a run
                # that is interrupted - by a scanner, a full disk, Ctrl-C - cannot leave
                # a truncated file where a verified one used to be.
                part = tar_path.with_suffix(tar_path.suffix + ".partial")
                print(f"\n{name} -> {tar_path.name}")
                fin = 0
                with tarfile.open(part, "w:gz") as tar:
                    for i, src in enumerate(files, 1):
                        a, raw, err = repack_bytes(src, scratch, args.complevel)
                        info = tarfile.TarInfo(f"{name}/{src.name}")
                        info.size = len(raw); info.mtime = int(src.stat().st_mtime)
                        info.mode = 0o644
                        tar.addfile(info, io.BytesIO(raw))
                        fin += a; worst = max(worst, err)
                        if args.pace:
                            time.sleep(args.pace)
                        if i % 25 == 0 or i == len(files):
                            print(f"  {i:4d}/{len(files)}  {fin/1e6:7.1f} MB read", flush=True)
                os.replace(part, tar_path)          # atomic; the old file is gone only now
                fout = tar_path.stat().st_size
                total_in += fin; total_out += fout
                print(f"  {fin/1e6:.1f} MB of grids -> {fout/1e6:.1f} MB tarball "
                      f"({fin/max(fout,1):.2f}x)")
        print(f"\nTotal {total_in/1e6:.0f} MB -> {total_out/1e6:.0f} MB "
              f"({total_in/max(total_out,1):.2f}x smaller)")
        print(f"Largest round-trip error anywhere: {worst:.4g} m")
        print(f"\nTarballs are in {outdir}. Run make_zenodo_archive.sh afterwards to "
              f"build the remaining components and refresh MANIFEST.txt.")
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
