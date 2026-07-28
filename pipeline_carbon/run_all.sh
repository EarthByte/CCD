#!/usr/bin/env bash
# ============================================================================
# run_all.sh — one command, non-interactive, for the whole carbonate + CO2
# chain, starting from the up-to-date AREA-WEIGHTED hybrid CCD. DM2026 only.
#
#   Step 0  regenerate area-weighted hybrid CCD + propagate to step 7 / paper [python]
#   Step 7  carbonate thickness grids (MIN / MEAN / MAX sed-rate)             [python]
#   Step 8  DM2026 carbonate volume/area stats + well-data ground-truth       [python]
#   Step 9  CO2 analysis notebooks (headless) + CO2 figures                   [jupyter]
#
# Nothing here runs in the cloud sandbox; run locally (pygplates, gplately,
# pygmt, jupyter/papermill) with the Alfonso_etal_2024_DM26 input grids present.
#
# Usage:
#   ./run_all.sh                 # everything, from CCD regeneration onward
#   ./run_all.sh --skip-ccd      # keep the current CCD; start at step 7
#   ./run_all.sh --skip-step7    # reuse existing grids; stats + CO2 only
#   ./run_all.sh --only-step9    # CO2 analysis + figures only
#   ./run_all.sh --no-caffeinate # do NOT auto-wrap in caffeinate (see below)
#
# Resume the CO2 notebooks after a fixed error (skip the ones that already ran):
#   CO2_START_FROM=02 ./run_all.sh --only-step9      # start at 02-Subducted-Carbon
# (name or prefix; see 03_run_co2_notebooks.sh). Combine with --only-step9 so
# steps 0/7/8 are not repeated.
#
# macOS overnight safety: on Darwin this script re-execs itself under
# `caffeinate -i` automatically, so the Mac won't idle-sleep mid-run and stall
# the (multi-hour) step 7 — a bare ./run_all.sh is always sleep-safe. Opt out
# with --no-caffeinate or NO_CAFFEINATE=1. On battery, clamshell/lid-close sleep
# is NOT overridden by caffeinate; keep the Mac on AC power for overnight runs.
# ============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# --- auto-caffeinate on macOS (re-exec once, then release on exit) ------------
for a in "$@"; do [ "$a" = "--no-caffeinate" ] && NO_CAFFEINATE=1; done
if [ "${NO_CAFFEINATE:-0}" != 1 ] && [ "${_CAFFEINATED:-0}" != 1 ] \
   && [ "$(uname -s)" = "Darwin" ] && command -v caffeinate >/dev/null 2>&1; then
  export _CAFFEINATED=1
  echo "[run_all] macOS: re-exec under 'caffeinate -i' (no idle sleep during run; --no-caffeinate to disable)."
  exec caffeinate -i "$0" "$@"
fi

source ./config.sh

DO0=1; DO7=1; DO8=1; DO9=1
for a in "$@"; do case "$a" in
  --skip-ccd)      DO0=0 ;;
  --skip-step7)    DO7=0 ;;
  --only-step9)    DO0=0; DO7=0; DO8=0 ;;
  --skip-step9)    DO9=0 ;;
  --no-caffeinate) ;;                     # handled above; no-op here
  *) echo "unknown arg: $a" >&2; exit 2 ;;
esac; done

if [ "$DO0" = 1 ]; then
  echo "########## STEP 0: area-weighted hybrid CCD -> downstream inputs ##########"
  "$PYTHON" 00_update_hybrid_ccd.py $UPDATE_CCD_ARGS
fi

if [ "$DO7" = 1 ]; then
  echo "########## STEP 7: carbonate thickness (min/mean/max) ##########"
  "$PYTHON" 01_carbonate_thickness_min_mean_max.py

  echo "---- handoff: link step-7 carbonate grids -> step-9 CarbonateSediment ----"
  # Step 9 (01-Sources-of-Carbon) reads decompacted grids from a NESTED layout:
  #   CarbonateSediment/carbonate_sed_thickness_DM2026/
  #     carbonate_sed_thickness_{min,mean,max}_DM2026/decompacted_sediment_thickness_0.25_<t>.nc
  # NB: step 7 writes the MEAN set to carbonate_sed_thickness_DM2026 (no _mean_),
  # so the mean is remapped to carbonate_sed_thickness_mean_DM2026 here.
  # Relative symlinks point back at the step-7 outputs: single source of truth,
  # no multi-GB duplication, and they survive the repo being relocated.
  CARB="$DM26_DIR/Grids/InputGrids/CarbonateSediment"
  PARENT="$CARB/carbonate_sed_thickness_DM2026"
  mkdir -p "$PARENT"
  # scrub any earlier FLAT staging so the nested layout is unambiguous
  find "$PARENT" -maxdepth 1 \( -name '*.nc' -o -name '*.xy' \) -delete 2>/dev/null || true
  rm -rf "$CARB/carbonate_sed_thickness_min_DM2026" "$CARB/carbonate_sed_thickness_max_DM2026"
  REL="$("$PYTHON" -c 'import os,sys; print(os.path.relpath(sys.argv[1], sys.argv[2]))' "$STEP7_DIR" "$PARENT")"
  link_scenario () {   # $1 = step-7 dir name, $2 = notebook subdir name
    if [ -d "$STEP7_DIR/$1" ]; then
      ln -sfn "$REL/$1" "$PARENT/$2"
      echo "  link $2 -> $REL/$1"
    else
      echo "  WARNING: step-7 output missing: $STEP7_DIR/$1" >&2
    fi
  }
  link_scenario carbonate_sed_thickness_min_DM2026  carbonate_sed_thickness_min_DM2026
  link_scenario carbonate_sed_thickness_DM2026      carbonate_sed_thickness_mean_DM2026
  link_scenario carbonate_sed_thickness_max_DM2026  carbonate_sed_thickness_max_DM2026
fi

if [ "$DO8" = 1 ]; then
  echo "########## STEP 8: DM2026 carbonate volume/area stats ##########"
  "$PYTHON" 02_carbonate_volume_stats.py --tmin "$STEP7_TIME_MIN" --tmax "$STEP7_TIME_MAX"
  echo "########## STEP 8: well-data ground-truth ##########"
  ( cd "$STEP8_DIR" && "$PYTHON" 08_ground_truth_carbonate_thickness.py )
fi

if [ "$DO9" = 1 ]; then
  echo "########## STEP 9: CO2 analysis notebooks (headless) ##########"
  bash ./03_run_co2_notebooks.sh
  echo "########## STEP 9: CO2 figures ##########"
  ( cd "$STEP9_DIR/standalone_plots" && "$PYTHON" run_all.py )
fi

echo
echo "Pipeline complete (from area-weighted CCD through CO2 analysis)."
