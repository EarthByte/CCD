#!/usr/bin/env bash
# ============================================================================
# 03_run_co2_notebooks.sh — headless execution of the step-9 CO2 notebooks.
#
# Each notebook is executed top-to-bottom, in order, with no interaction, using
# `jupyter nbconvert --to notebook --execute --inplace` (equivalently papermill).
# The notebook file is updated in place with the executed outputs, and any CSVs /
# grids / figures the cells write are produced exactly as in an interactive run.
#
# Requires: jupyter + nbconvert (or papermill) in the active environment, plus
# all step-9 dependencies (gplately, pygplates, pygmt, ...). Cannot run in the
# cloud sandbox.
# ============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

run_nb () {
  local nb="$1"
  local path="$CO2_NB_DIR/$nb"
  if [ ! -f "$path" ]; then echo "  SKIP (missing): $path"; return 0; fi
  echo "==> executing: $nb"
  if command -v papermill >/dev/null 2>&1; then
    papermill "$path" "$path" -k "$JUPYTER_KERNEL" --no-progress-bar
  else
    "$PYTHON" -m jupyter nbconvert --to notebook --execute --inplace \
      --ExecutePreprocessor.timeout="$NB_TIMEOUT" \
      --ExecutePreprocessor.kernel_name="$JUPYTER_KERNEL" "$path"
  fi
}

# Optional resume point: notebook name or prefix (env CO2_START_FROM, or $1).
# Notebooks before it are skipped so you can continue after a fixed error without
# re-running the ones that already succeeded. Example:
#   CO2_START_FROM=02 bash 03_run_co2_notebooks.sh
#   bash 03_run_co2_notebooks.sh 02-Subducted-Carbon.ipynb
START="${CO2_START_FROM:-${1:-}}"
started=1
[ -n "$START" ] && started=0

run_list () {
  for nb in "$@"; do
    [ -z "$nb" ] && continue
    if [ "$started" = 0 ]; then
      case "$nb" in
        "$START"*) started=1 ;;
        *) echo "  skip (resume: before $START): $nb"; continue ;;
      esac
    fi
    run_nb "$nb"
  done
}

echo "### Step 9: CO2 analysis (headless notebooks) in $CO2_NB_DIR"
[ -n "$START" ] && echo "### resuming from: $START"
run_list "${CO2_PREP_NOTEBOOKS[@]:-}"
run_list "${CO2_NOTEBOOKS[@]}"
if [ "$started" = 0 ]; then
  echo "ERROR: resume notebook '$START' not found in the notebook lists" >&2
  exit 3
fi
echo "### Step 9 notebooks complete."
