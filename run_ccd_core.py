#!/usr/bin/env python3
"""
Runner for the CCD core — steps 1-7 (regional CCD -> global synthesis ->
sea-level calibration -> planktogenic hybrid CCD). Run from the repo root:

    python run_ccd_core.py

Each step reads the declared outputs of earlier steps (see ccdworkflow/config.py)
and writes into its own steps/<step>/outputs/ folder. Step 2 (ocean-basin areas)
is a notebook prep that generates the basin-area fractions consumed by step 3 and
is not part of this automated numeric core. The carbonate + carbon-cycle stages
(steps 8-10) are run separately with ./run_carbon_cycle.sh.
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

STEPS = [
    ("steps/step1_regional_ccd_resampling/resample_regional_ccd.py", []),
    ("steps/step3_global_ccd_synthesis/global_ccd_synthesis.py", []),
    ("steps/step4_sealevel_envelope/sealevel_envelope.py", ["--validate"]),
    ("steps/step5_ccd_lowpass_filter/ccd_lowpass_filter.py", []),
    ("steps/step6_sealevel_ccd_regression/regression_and_prediction.py", []),
    ("steps/step7_planktogenic_hybrid_ccd/planktogenic_hybrid.py", []),
]


def main() -> None:
    root = Path(__file__).resolve().parent
    sys.path.insert(0, str(root))
    for s, extra in STEPS:
        print(f"\n===== {s} =====")
        sys.argv = [str(root / s)] + extra
        runpy.run_path(str(root / s), run_name="__main__")
    print("\nCCD core complete (steps 1-7). Outputs in each steps/<step>/outputs/.")


if __name__ == "__main__":
    main()
