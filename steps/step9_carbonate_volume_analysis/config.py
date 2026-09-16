"""
Shared configuration for the carbonate-thickness comparison analysis.

Two grid sets:
    DM2026 -- compacted carbonate sediment thickness using the Dutkiewicz &
              Muller (2026, in prep) global CCD curve.
    BW1991 -- compacted carbonate sediment thickness using the Boss &
              Wilkinson (1991) CCD curve.

Both products use the same paleobathymetry, computed on the Alfonso et al. (2024)
plate model following the pyBacktrack 1.5 method of Müller et al. (2026),
as input.  Grids are 0.25 deg lat/lon, present-day grid (no plate
reconstruction applied for plotting), one NetCDF per Ma from 0 to 170 Ma.
Variable name in each file: ``z`` (metres of compacted carbonate sediment).
NaN cells = non-carbonate-bearing (above CCD, on land, or outside the
deposition mask).

This module is imported by every downstream script in this folder.
"""
from __future__ import annotations

import functools
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import xarray as xr


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent.parent      # .../Dutkiewicz_Muller_CCD
GRID_ROOT = HERE.parent / "step7_carbonate"

DM2026_DIR = GRID_ROOT / "carbonate_sed_thickness_DM2026"
BW1991_DIR = GRID_ROOT / "carbonate_sed_thickness_BW1991"

# Alfonso et al. 2024 / pyBacktrack merged paleobathymetry input grids.
# Used both as the input to the carbonate thickness workflow and (here) as
# a standalone video product. One NetCDF per integer Ma from 0 to 170 Ma.
PALEOBATH_DIR = GRID_ROOT / "input_grids" / "Alfonso2024_pybacktrack_merged_paleobathymetry"
PALEOBATH_FMT = "paleobathymetry_{t:d}Ma.nc"

OUTPUT_DIR = HERE / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
(OUTPUT_DIR / "stats").mkdir(exist_ok=True)
(OUTPUT_DIR / "figures").mkdir(exist_ok=True)
(OUTPUT_DIR / "video_frames").mkdir(exist_ok=True)

# Filename template: <type>_0.25_<time>.nc
COMPACTED_FMT = "compacted_sediment_thickness_0.25_{t:d}.nc"

# Variable name inside each NetCDF
GRID_VAR = "z"

# ---------------------------------------------------------------------------
# Time range — user asked for 0..150 Ma (carbonate was negligible before 150 Ma).
# ---------------------------------------------------------------------------
MIN_TIME_MA = 0
MAX_TIME_MA = 150
TIME_STEP_MYR = 1
TIMES = list(range(MIN_TIME_MA, MAX_TIME_MA + 1, TIME_STEP_MYR))   # 0..150 inclusive

# ---------------------------------------------------------------------------
# Boxplot binning (mirrors Fig 11 of the pyBacktrack repo)
# ---------------------------------------------------------------------------
BIN_WIDTH_MYR = 5
BIN_EDGES = np.arange(MIN_TIME_MA, MAX_TIME_MA + BIN_WIDTH_MYR, BIN_WIDTH_MYR)   # 0,5,...,150
BIN_CENTRES = (BIN_EDGES[:-1] + BIN_EDGES[1:]) / 2.0                              # 2.5,7.5,...,147.5

# ---------------------------------------------------------------------------
# Greedy non-maximum-suppression for the picked-times analysis.
# ---------------------------------------------------------------------------
N_PICKED_TIMES = 4
NMS_EXCLUSION_MYR = 25

# ---------------------------------------------------------------------------
# Latitude bands (10 deg) used by the region-characterisation analysis.
# ---------------------------------------------------------------------------
LAT_BAND_EDGES = np.arange(-90, 91, 10)              # -90, -80, ..., 90
LAT_BAND_CENTRES = (LAT_BAND_EDGES[:-1] + LAT_BAND_EDGES[1:]) / 2.0


# ---------------------------------------------------------------------------
# Present-day ocean-basin classifier.
#
# Simple lon/lat cuts — coarse but transparent and reproducible without
# requiring shapefile data.  Order matters: a cell is assigned to the first
# basin whose predicate returns True.  Adjust the cuts here if you want
# finer-grained polygons later.
# ---------------------------------------------------------------------------
BASINS = ("Arctic", "Southern", "Pacific", "Atlantic", "Indian")


def basin_of(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    """Return a string-dtype array of basin names for each (lat, lon) cell.

    `lat_deg` and `lon_deg` must be broadcastable to a common shape.

    Cuts (matched in this order so polar / circum-Antarctic dominate the
    high-latitude assignments before the basin lon-cuts kick in):

        Arctic            : lat >= 66
        Southern          : lat <= -60
        Pacific           : lon < -70  OR  lon >= 145
        Atlantic          : -70 <= lon < 20
        Indian            : 20 <= lon < 145

    lon is assumed to be in [-180, 180].
    """
    lat = np.broadcast_to(lat_deg, np.broadcast_shapes(lat_deg.shape, lon_deg.shape))
    lon = np.broadcast_to(lon_deg, lat.shape)
    out = np.empty(lat.shape, dtype=object)
    out[:] = ""
    arctic = lat >= 66
    southern = (lat <= -60) & ~arctic
    pacific = ~(arctic | southern) & ((lon < -70) | (lon >= 145))
    atlantic = ~(arctic | southern | pacific) & (lon >= -70) & (lon < 20)
    indian = ~(arctic | southern | pacific | atlantic) & (lon >= 20) & (lon < 145)
    out[arctic] = "Arctic"
    out[southern] = "Southern"
    out[pacific] = "Pacific"
    out[atlantic] = "Atlantic"
    out[indian] = "Indian"
    return out


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------
def grid_path(source: str, t: int) -> Path:
    """source in {'DM2026', 'BW1991'}; t in 0..170 Ma."""
    if source == "DM2026":
        return DM2026_DIR / COMPACTED_FMT.format(t=t)
    if source == "BW1991":
        return BW1991_DIR / COMPACTED_FMT.format(t=t)
    raise ValueError(f"Unknown source {source!r}")


@functools.lru_cache(maxsize=64)
def load_grid(source: str, t: int) -> xr.DataArray:
    """Load one (source, time) compacted-carbonate thickness grid as an
    xarray DataArray with coords (lat, lon) and NaN-masked non-deposition
    cells.  Cached so the boxplot / stats / video scripts don't pay the
    open-NetCDF cost twice per call site.
    """
    path = grid_path(source, t)
    if not path.exists():
        raise FileNotFoundError(path)
    return xr.open_dataset(path)[GRID_VAR].load()


@functools.lru_cache(maxsize=64)
def load_paleobath(t: int) -> xr.DataArray:
    """Load the Alfonso 2024 / pyBacktrack merged paleobathymetry grid at
    age ``t`` Ma as an xarray DataArray. Values are depth in metres
    (negative below sea level; NaN over land / un-reconstructed cells).
    Cached as for load_grid().
    """
    path = PALEOBATH_DIR / PALEOBATH_FMT.format(t=t)
    if not path.exists():
        raise FileNotFoundError(path)
    ds = xr.open_dataset(path)
    # The NetCDF can carry the data variable under different names depending
    # on the producer ("z", "elevation", "paleobathymetry", ...). Pick the
    # first data_var rather than hard-coding so this stays robust.
    var = list(ds.data_vars)[0]
    return ds[var].load()


def common_mask(t: int) -> np.ndarray:
    """Boolean mask of cells where BOTH products have finite carbonate
    thickness at age t.  Use this for the difference grid + per-time
    statistics so we only ever compare apples to apples.
    """
    a = load_grid("DM2026", t).values
    b = load_grid("BW1991", t).values
    return np.isfinite(a) & np.isfinite(b)


# ---------------------------------------------------------------------------
# Registration normaliser
#
# All grids in this workflow are produced and consumed as *gridline*-
# registered global grids on the -180..180 longitude convention (matching
# what `carbonate_sediment_thickness.py` writes via `gmt nearneighbor -rg`
# in -R-180/180/-90/90). If any input grid sneaks through pixel-registered,
# GMT's grdimage emits:
#
#   "Longitude range too small; geographic boundary condition changed
#    to natural."
#
# ...and falls back from periodic to natural boundary conditions, which can
# subtly distort the western/eastern edges of the rendered map. This helper
# materialises an xarray.DataArray to a temporary NetCDF and runs
# `gmt grdedit` to *guarantee* gridline registration before the grid is
# handed to pygmt's grdimage. Cheap (a single grdedit per frame) and
# completely silences the warning regardless of how the grid was produced.
# ---------------------------------------------------------------------------
def to_gridline_nc(da: xr.DataArray) -> Path:
    """Write `da` to a temp NetCDF and force gridline registration.

    Returns the temp file path. Caller is responsible for unlinking it
    (or letting the OS clean the temp dir).
    """
    tmp = Path(tempfile.mkstemp(suffix=".nc", prefix="ccd_gridline_")[1])
    # xarray writes a CF-compliant NetCDF; GMT will pick up the spatial coords.
    da.to_netcdf(tmp)
    # Probe the current registration. `gmt grdinfo` prints e.g.
    #   "Gridline node registration used"   or
    #   "Pixel node registration used"
    info = subprocess.run(
        ["gmt", "grdinfo", str(tmp)],
        capture_output=True, text=True, check=True
    ).stdout
    if "Pixel node registration used" in info:
        # `grdedit -T` toggles registration in place (pixel -> gridline here).
        subprocess.run(
            ["gmt", "grdedit", str(tmp), "-T"],
            check=True
        )
    return tmp


def diff_grid(t: int) -> xr.DataArray:
    """DM2026 - BW1991 at age t, NaN wherever EITHER input is NaN.

    Built with explicit numpy ops (not xarray.where) so the NaN footprint
    is bulletproof against any dtype / mask-propagation surprises in
    downstream pyGMT calls (a previous version used xarray.where and the
    delta video painted continents with the cpt background colour).

    Positive = DM2026 produces more carbonate than BW1991.
    """
    a_da = load_grid("DM2026", t)
    b_da = load_grid("BW1991", t)
    a = a_da.values
    b = b_da.values
    d = (a - b).astype(np.float32)
    bad = ~(np.isfinite(a) & np.isfinite(b))
    d[bad] = np.float32("nan")
    return xr.DataArray(d,
                        coords={"lat": a_da["lat"].values,
                                "lon": a_da["lon"].values},
                        dims=("lat", "lon"),
                        name="delta")


# ---------------------------------------------------------------------------
# Print banner on import so every figure / video run logs the active config.
# Opt out by setting CCD_QUIET=1 in the environment.
# ---------------------------------------------------------------------------
if not os.environ.get("CCD_QUIET"):
    print(f"[ccd-cmp] DM2026 dir : {DM2026_DIR}")
    print(f"[ccd-cmp] BW1991 dir : {BW1991_DIR}")
    print(f"[ccd-cmp] times      : {MIN_TIME_MA}..{MAX_TIME_MA} Ma "
          f"in {TIME_STEP_MYR} Myr steps ({len(TIMES)} frames)")
    print(f"[ccd-cmp] output dir : {OUTPUT_DIR}")
