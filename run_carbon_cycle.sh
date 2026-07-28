#!/usr/bin/env bash
# Top-level runner for the carbonate + carbon-cycle stages (steps 8-10):
#   step 8  carbonate sediment thickness (min/mean/max)
#   step 9  carbonate volume/area stats + well-data ground-truth
#   step 10 CO2 / carbon-cycle notebooks (headless) + figures
# Thin wrapper around pipeline_carbon/run_all.sh (which also regenerates the
# hybrid CCD from run_ccd_core.py first). Runs locally, not in the cloud.
exec "$(dirname "${BASH_SOURCE[0]}")/pipeline_carbon/run_all.sh" "$@"
