"""
Central configuration: directory layout and canonical parameters.

All paths derive from the repository root so the workflow is fully relocatable.
Import ``REPO_ROOT``, ``DATA``, the per-step output dirs and the canonical file
constants rather than hard-coding paths inside step scripts.

Layout (steps 1-10):
    steps/step1_timescale_conversion/        -> STEP1_DIR (outputs/); every CCD and
                                              sea-level curve normalised to GTS2020
    steps/step2_ocean_basin_areas/           (notebook prep; produces the area
                                              fractions used by step 3)
    steps/step3_global_ccd_synthesis/        -> STEP3_DIR
    steps/step4_sealevel_envelope/           -> STEP4_DIR
    steps/step5_ccd_lowpass_filter/          -> STEP5_DIR
    steps/step6_sealevel_ccd_regression/     -> STEP6_DIR
    steps/step7_planktogenic_hybrid_ccd/     -> STEP7_DIR
    steps/step8_carbonate_sediment_thickness/  (carbon runner)
    steps/step9_carbonate_volume_analysis/     (carbon runner)
    steps/step10_carbon_cycle_degassing/       (carbon runner)
Each step folder holds its own ``outputs/`` (csv + figures/).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------
# Directory layout
# --------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent

DATA = REPO_ROOT / "data"                       # raw inputs (not produced by any step)
REGIONAL_CCD = DATA / "regional_ccd"            # archived 0.5 Myr ATL/IND/PAC curves;
                                                # kept for provenance, read by no step
SEALEVEL = DATA / "sealevel"                    # sea-level compilations
REFERENCE_CCD = DATA / "reference_ccd"          # published reference CCD curves
CCD_FRACTIONS_XLSX = DATA / "ccds_and_oceanbasin_fractions.xlsx"

STEPS = REPO_ROOT / "steps"
FIGURES = REPO_ROOT / "Figures"                 # collected diagnostic figure set (PNG + PDF)


def _paper_figures_dir() -> Path | None:
    """Folder holding the paper's numbered figures, or None if it cannot be found.

    Figure 1 is the only paper figure produced by a pipeline step rather than by a
    make_fig*.py script, so step 3 writes a copy straight into it. Two layouts are
    supported: the working tree (this repo sits next to Paper/) and the public repo
    (this repo IS the repo root and the figures live in figures/). Override with
    the CCD_PAPER_FIGURES environment variable.
    """
    env = os.environ.get("CCD_PAPER_FIGURES")
    if env:
        return Path(env).expanduser()
    for cand in (REPO_ROOT.parent / "Paper" / "Figures", REPO_ROOT / "figures"):
        if cand.is_dir():
            return cand
    return None


PAPER_FIGURES = _paper_figures_dir()

# Per-step output directories (each inside its step folder: steps/<step>/outputs)
STEP1_DIR = STEPS / "step1_timescale_conversion" / "outputs"        # GTS2020-normalised inputs
STEP3_DIR = STEPS / "step3_global_ccd_synthesis" / "outputs"        # synthesised global CCD
STEP4_DIR = STEPS / "step4_sealevel_envelope" / "outputs"           # sea-level quantile envelope
STEP5_DIR = STEPS / "step5_ccd_lowpass_filter" / "outputs"          # low-pass filtered CCD
STEP6_DIR = STEPS / "step6_sealevel_ccd_regression" / "outputs"     # regression + prediction
STEP7_DIR = STEPS / "step7_planktogenic_hybrid_ccd" / "outputs"     # planktogenic hybrid CCD

ALL_STEP_DIRS = [STEP1_DIR, STEP3_DIR, STEP4_DIR, STEP5_DIR, STEP6_DIR, STEP7_DIR]

# Short step key -> its output directory (used by the figure savers to put each
# step's figures in its own outputs/figures/ folder, in addition to Figures/).
STEP_DIRS = {
    "step1": STEP1_DIR,
    "step3": STEP3_DIR,
    "step4": STEP4_DIR,
    "step5": STEP5_DIR,
    "step6": STEP6_DIR,
    "step7": STEP7_DIR,
}


def step_figures_dir(step: str) -> Path:
    """Per-step figures folder: ``steps/<step>/outputs/figures``."""
    base = STEP_DIRS.get(step, STEPS / step / "outputs")
    return base / "figures"


def ensure_dirs() -> None:
    """Create the per-step output dirs and the Figures folder if missing."""
    for d in ALL_STEP_DIRS:
        d.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Canonical output file names (single source of truth for step-to-step wiring)
# --------------------------------------------------------------------------
# Step 1 -> Steps 3 and 4  (everything normalised to GTS2020)
GTS2020_ATL_CCD = STEP1_DIR / "ATL_CCD_GTS2020_1my.txt"
GTS2020_PAC_CCD = STEP1_DIR / "PAC_CCD_GTS2020_1my.txt"
GTS2020_IND_CCD = STEP1_DIR / "IND_CCD_GTS2020_1my.txt"
GTS2020_REGIONAL_CCD = {"Atlantic": GTS2020_ATL_CCD,
                        "Pacific": GTS2020_PAC_CCD,
                        "Indian": GTS2020_IND_CCD}
GTS2020_SL_SHORTTERM = STEP1_DIR / "sealevel_shortterm_hybrid_GTS2020.txt"
GTS2020_HAQ_SHORTTERM = STEP1_DIR / "Haq_shortterm_hybrid_GTS2020.txt"
GTS2020_HAQ_LONGTERM = STEP1_DIR / "Haq87_longterm_GTS2020.txt"

# Step 3 -> Step 5
GLOBAL_CCD_0_52 = STEP3_DIR / "global_ccd_with_basin_dispersion_0-52Ma.txt"

# Step 4 -> Steps 5, 6
SL_ENVELOPE_FULL = STEP4_DIR / "sea_level_quantile_envelope_0-205Ma.txt"   # full range
SL_ENVELOPE_0_52 = STEP4_DIR / "sea_level_quantile_envelope_0-52Ma.txt"    # calibration subset

# Step 5 -> Steps 6, 7
LOWPASS_CCD_0_52 = STEP5_DIR / "low_pass_filtered_ccd_0-52Ma.txt"

# Step 6 -> Step 7
PREDICTED_CCD_SL = STEP6_DIR / "predicted_ccd_sl_0-205Ma.txt"

# Step 7 (final products)
HYBRID_CCD = STEP7_DIR / "hybrid_ccd_obs_pred_combined.txt"
HYBRID_CCD_ENVELOPE = STEP7_DIR / "hybrid_ccd_error_envelope.txt"


# --------------------------------------------------------------------------
# Canonical numerical parameters
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    """Fixed parameters for the pipeline (pinned so runs are deterministic)."""

    # ---- sea-level rolling high-quantile envelope ---------------------------
    sl_window_myr: float = 2.0
    sl_quantile: float = 0.90
    sl_smooth_myr: float = 1.5
    sl_grid_dt: float = 0.1

    # ---- CCD Butterworth low-pass (matched to SL bandwidth) ----------------
    ccd_butter_cutoff: float = 0.5
    ccd_butter_order: int = 1
    ccd_smooth_samples: int = 1
    ccd_fs: float = 1.0
    ccd_cal_max_age: float = 52.0

    # ---- preferred regression (reduced major axis) -------------------------
    use_fitted_rma: bool = True
    pred_slope: float = 5.72
    pred_intercept: float = -4456.8
    pred_tstart: float = 0.0
    pred_tend: float = 205.0
    pred_dt: float = 1.0
    pred_corr_lengthscale_myr: float = 10.0
    pred_z95: float = 1.96

    # ---- planktogenic correction + splice ----------------------------------
    plankton_slope_m_per_myr: float = 20.0
    plankton_start_age_ma: float = 115.0
    splice_obs_max_age: float = 52.0
    splice_bridge_age: float = 53.0
    splice_pred_min_age: float = 54.0


PARAMS = Params()
