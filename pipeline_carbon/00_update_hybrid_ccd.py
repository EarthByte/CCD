#!/usr/bin/env python3
"""
Step 0 (batch) — regenerate the up-to-date AREA-WEIGHTED hybrid CCD and push it
into the downstream inputs, so nothing stale propagates through steps 7–9.

What it does:
  1. Runs the clean CCD workflow (CCD_workflow_clean/run_ccd_core.py, steps 1–7) which
     produces the area-weighted hybrid CCD.
  2. Reformats that curve to the 0–170 Ma, 4-column, space-separated form that the
     carbonate-thickness code expects (age  mean  min  max; CCD negative-down).
  3. Writes it to BOTH consumers, overwriting the older equal-weight versions:
       steps/step8_carbonate_sediment_thickness/input_data/CCD_sl_hybrid_2026.txt  (step 8)
       Paper/Figures/CCD_hybrid_DM2026.txt                              (paper/attribution)

Run this first (run_all.sh does) so step 7 rebuilds carbonate grids from the
current CCD, step 8 restats them, and step 9 re-derives subducted carbon — the
whole chain stays consistent with the latest area-weighted CCD.

The clean workflow (numpy/pandas/scipy/matplotlib) runs anywhere; steps 7–9 need
your geoscience environment.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CLEAN_REPO = Path(os.environ.get("CLEAN_REPO", HERE.parent))    # CCD_workflow_clean
PROJECT = Path(os.environ.get("PROJECT_ROOT", CLEAN_REPO.parent))  # Dutkiewicz_Muller_CCD
STEPS_ROOT = Path(os.environ.get("STEPS_ROOT", CLEAN_REPO / "steps"))
# Step 7 writes into its own steps/ folder; the old outputs/step6_hybrid/ path is
# from the pre-renumbering layout and no longer exists, which silently left
# Paper/Figures/CCD_hybrid_DM2026.txt stale. Fall back to the old path only if
# someone is running an archived tree.
HYBRID_SRC = STEPS_ROOT / "step7_planktogenic_hybrid_ccd" / "outputs" / "hybrid_ccd_obs_pred_combined.txt"
if not HYBRID_SRC.exists():
    _legacy = CLEAN_REPO / "outputs" / "step6_hybrid" / "hybrid_ccd_obs_pred_combined.txt"
    if _legacy.exists():
        HYBRID_SRC = _legacy

CONSUMERS = [
    STEPS_ROOT / "step8_carbonate_sediment_thickness" / "input_data" / "CCD_sl_hybrid_2026.txt",
    PROJECT / "Paper" / "Figures" / "CCD_hybrid_DM2026.txt",
]


def regenerate_clean_ccd() -> None:
    if not (CLEAN_REPO / "run_ccd_core.py").exists():
        raise SystemExit(f"clean workflow not found at {CLEAN_REPO} "
                         f"(set CLEAN_REPO env var to its path)")
    print(f"[step0] running clean CCD workflow: {CLEAN_REPO}/run_ccd_core.py")
    subprocess.run([sys.executable, "run_ccd_core.py"], cwd=str(CLEAN_REPO), check=True,
                   env={**os.environ, "PYTHONPATH": str(CLEAN_REPO)})


def reformat_to_step7(src: Path) -> str:
    df = pd.read_csv(src, sep=r"\s+|\t+", engine="python", comment="#", header=None,
                     names=["age", "mean", "mn", "mx"])
    df = df[(df["age"] >= 0) & (df["age"] <= 170)]
    grid = np.arange(0, 171)
    out = {c: np.interp(grid, df["age"], df[c]) for c in ["mean", "mn", "mx"]}
    lines = [f"{int(a)} {out['mean'][i]:.2f} {out['mn'][i]:.2f} {out['mx'][i]:.2f}"
             for i, a in enumerate(grid)]
    return "\n".join(lines) + "\n"


def main() -> None:
    if "--no-regen" not in sys.argv:
        regenerate_clean_ccd()
    if not HYBRID_SRC.exists():
        raise SystemExit(f"hybrid CCD not found: {HYBRID_SRC}")
    text = reformat_to_step7(HYBRID_SRC)
    for dest in CONSUMERS:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(text)
        print(f"[step0] updated {dest}")
    print("[step0] area-weighted hybrid CCD propagated to step 7 + paper figures.")


if __name__ == "__main__":
    main()
