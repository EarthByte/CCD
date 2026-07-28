#!/usr/bin/env python3
"""
Step 7 (batch) — carbonate sediment thickness grids for MIN / MEAN / MAX.

This reproduces ``carbonate_sediment_thickness_2026.ipynb`` non-interactively and
runs all three scenarios in one go. In that notebook the three scenarios are
selected by the **maximum carbonate sedimentation-rate curve** (not by the CCD);
the CCD is fixed at the DM2026 hybrid curve, and the model choice (DM2026 vs
BW1991) is the CCD curve. So:

    scenario   sed-rate curve            output directory
    --------   -----------------------   --------------------------------
    min        sed_rate_min.txt          carbonate_sed_thickness_min_DM2026
    mean       sed_rate_best.txt         carbonate_sed_thickness_DM2026
    max        sed_rate_max.txt          carbonate_sed_thickness_max_DM2026

with the CCD fixed at input_data/CCD_sl_hybrid_2026.txt for all three.

The driver chdir's into the step-7 directory so all relative paths (input_data/,
input_grids/, plate-model-repo/, output dirs) resolve exactly as in the notebook.

CANNOT run in the cloud sandbox: needs pygplates, GMT, the Alfonso2024 plate
model, and the DM26 age + paleobathymetry grids. Run in your local environment.
Parameters below mirror the notebook one-to-one — change them there and here
together if you retune.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_STEPS_ROOT = Path(os.environ.get("STEPS_ROOT",
                   Path(__file__).resolve().parent.parent / "steps"))
STEP7 = _STEPS_ROOT / "step8_carbonate_sediment_thickness"
os.chdir(STEP7)                      # match the notebook's working directory
sys.path.insert(0, str(STEP7))
import carbonate_sediment_thickness   # noqa: E402

# ---------------------------------------------------------------------------
# Fixed parameters (identical to carbonate_sediment_thickness_2026.ipynb)
# ---------------------------------------------------------------------------
GRID_SPACING = 0.25
MIN_LAT, MAX_LAT = -90, 90
MIN_LON, MAX_LON = -180, 180
TIMES = range(0, 171)
USE_ALL_CPU_CORES = True
CARBONATE_ANCHOR_PLATE_ID = 0

# Plate model: Alfonso2024 via plate-model-manager (downloads to plate-model-repo),
# matching the notebook (using_local_model = False). Set USE_LOCAL_MODEL = True to
# use a bundled local topology model instead (offline).
USE_LOCAL_MODEL = False
import glob
if USE_LOCAL_MODEL:
    ROTATION_FILENAMES = glob.glob(os.path.join("input_data", "topology_model", "2019_v2", "*.rot"))
    TOPOLOGY_FILENAMES = glob.glob(os.path.join("input_data", "topology_model", "2019_v2", "*.gpmlz"))
else:
    from plate_model_manager import PlateModelManager
    _pm = PlateModelManager().get_model("Alfonso2024", data_dir="plate-model-repo")
    ROTATION_FILENAMES = _pm.get_rotation_model()
    TOPOLOGY_FILENAMES = _pm.get_topologies()

# Fixed DM2026 hybrid CCD curve (time -> CCD, negative).
CCD_CURVE = "input_data/CCD_sl_hybrid_2026.txt"

# Age + bathymetry grids (Alfonso 2024, 2026 build).
AGE_GRID_FORMAT = "./input_grids/Alfonso2024_SeafloorAgeGrids-2026/SEAFLOOR_AGE_grid_{:.2f}Ma.nc"
AGE_GRID_ANCHOR_PLATE_ID = 0
BATHYMETRY_GRID_FORMAT = "./input_grids/Alfonso2024_pybacktrack_merged_paleobathymetry/paleobathymetry_{:.0f}Ma.nc"
BATHYMETRY_GRID_ANCHOR_PLATE_ID = 0
BATHYMETRY_OLDEST_TIME = 170

# scenario -> (max carbonate sed-rate curve, output directory)
SCENARIOS = {
    "min":  ("input_data/sed_rate_min.txt",  "carbonate_sed_thickness_min_DM2026"),
    "mean": ("input_data/sed_rate_best.txt", "carbonate_sed_thickness_DM2026"),
    "max":  ("input_data/sed_rate_max.txt",  "carbonate_sed_thickness_max_DM2026"),
}


def run_scenario(name: str, sed_rate_curve: str, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n===== step 7: carbonate thickness [{name}]  sed_rate={sed_rate_curve} -> {out_dir} =====")
    carbonate_sediment_thickness.calc_sedimentation_and_write_data_for_times(
        TIMES,
        (MIN_LAT, MAX_LAT),
        (MIN_LON, MAX_LON),
        GRID_SPACING,
        ROTATION_FILENAMES,
        TOPOLOGY_FILENAMES,
        CCD_CURVE,
        sed_rate_curve,
        AGE_GRID_FORMAT,
        AGE_GRID_ANCHOR_PLATE_ID,
        BATHYMETRY_GRID_FORMAT,
        BATHYMETRY_GRID_ANCHOR_PLATE_ID,
        BATHYMETRY_OLDEST_TIME,
        os.path.join(out_dir, "decompacted_sediment_thickness"),
        os.path.join(out_dir, "compacted_sediment_thickness"),
        os.path.join(out_dir, "deposition_mask"),
        CARBONATE_ANCHOR_PLATE_ID,
        USE_ALL_CPU_CORES,
    )


def main() -> None:
    only = sys.argv[1:] or ["min", "mean", "max"]     # e.g. `python 01_...py mean`
    for name in only:
        sed_rate_curve, out_dir = SCENARIOS[name]
        run_scenario(name, sed_rate_curve, out_dir)
    print("\nStep 7 complete: min/mean/max carbonate thickness grids written.")


if __name__ == "__main__":
    main()
