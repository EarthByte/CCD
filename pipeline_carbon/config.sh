#!/usr/bin/env bash
# ============================================================================
# config.sh — configuration for the carbonate + CO2 batch pipeline (steps 7-9).
# Lives inside CCD_workflow_clean/pipeline_carbon/. DM2026 only.
#
# EDIT the paths marked  # >>> CONFIRM.  Nothing here runs in the cloud sandbox;
# run it locally with pygplates / gplately / GMT / jupyter and the input grids.
# ============================================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLEAN_REPO="$(cd "$PIPELINE_DIR/.." && pwd)"                # CCD_workflow_clean
PROJECT_ROOT="$(cd "$CLEAN_REPO/.." && pwd)"               # Dutkiewicz_Muller_CCD

# Where the step 7/8/9 code + input grids live. Consolidated into this repo:
# defaults to $CLEAN_REPO/steps (step8_carbonate_sediment_thickness,
# step9_carbonate_volume_analysis, step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26).
# Override only if you relocate the heavy grids elsewhere.
export STEPS_ROOT="${STEPS_ROOT:-$CLEAN_REPO/steps}"   # in-repo (consolidated, steps 8-10)
export CLEAN_REPO PROJECT_ROOT

STEP7_DIR="$STEPS_ROOT/step8_carbonate_sediment_thickness"
STEP8_DIR="$STEPS_ROOT/step9_carbonate_volume_analysis"
STEP9_DIR="$STEPS_ROOT/step10_carbon_cycle_degassing"
DM26_DIR="$STEP9_DIR/Alfonso_etal_2024_DM26"

# Python for your geoscience env (pygplates, gplately, pygmt, xarray, jupyter).
PYTHON="${PYTHON:-python}"

# ---- Step 0: regenerate + propagate the area-weighted hybrid CCD ------------
# 00_update_hybrid_ccd.py runs CCD_workflow_clean/run_ccd_core.py and writes the
# area-weighted CCD into step 7 input + Paper/Figures. Set --no-regen to skip
# re-running the clean workflow and just re-propagate the existing curve.
UPDATE_CCD_ARGS="${UPDATE_CCD_ARGS:-}"

# ---- Step 7: carbonate sediment thickness (min/mean/max = sed-rate curve) ----
STEP7_OUT_MEAN="$STEP7_DIR/carbonate_sed_thickness_DM2026"
STEP7_OUT_MIN="$STEP7_DIR/carbonate_sed_thickness_min_DM2026"
STEP7_OUT_MAX="$STEP7_DIR/carbonate_sed_thickness_max_DM2026"
STEP7_TIME_MIN=0
STEP7_TIME_MAX=170

# ---- Step 9: CO2 notebooks (headless execution order) -----------------------
CO2_NOTEBOOKS=(
  "01-Sources-of-Carbon.ipynb"
  "02-Subducted-Carbon.ipynb"
  "03-Carbon-Degassing.ipynb"
  "04-Carbonate-Platform-Degassing.ipynb"
  "05-Atmospheric-Carbon.ipynb"
  "06-Plate-Tectonic-Stats.ipynb"
)
CO2_PREP_NOTEBOOKS=(
  # "utils/0E-MakeAgeSRGrids.ipynb"
  # "utils/Min-Mean-Max-Crustal-Carbon.ipynb"
)
CO2_NB_DIR="$DM26_DIR"
JUPYTER_KERNEL="${JUPYTER_KERNEL:-python3}"
NB_TIMEOUT="${NB_TIMEOUT:--1}"

# Resume the CO2 notebook sequence from a given notebook (name or prefix, e.g.
# "02" or "02-Subducted-Carbon.ipynb"); notebooks before it are skipped. Empty =
# run the whole list. Exported so `bash 03_run_co2_notebooks.sh` (a child of
# run_all.sh) inherits it. Override on the command line, e.g.:
#   CO2_START_FROM=02 ./run_all.sh --only-step9
export CO2_START_FROM="${CO2_START_FROM:-}"
