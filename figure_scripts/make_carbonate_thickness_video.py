#!/usr/bin/env python3
"""
Supplementary Movie S? -- compacted carbonate sediment thickness, 0-170 Ma.

Mollweide projection, coloured with the SAME palette as Figure 3
(data/carbonate_thickness_roma_pale.cpt), with GPlately-reconstructed
plate boundaries (black + white halo), subduction teeth, gray continental
polygons + thin coastlines. The carbonate-thickness grids are already in
Alfonso et al. (2025) palaeo-coordinates at each Ma, so they are plotted
directly (no raster reconstruction); only the vector overlays are
reconstructed to each frame's age.

Usage:
    python3 make_carbonate_thickness_video.py                 # 0-170 Ma, 1 Myr, 8 fps
    python3 make_carbonate_thickness_video.py --time-step 5   # quick preview
    python3 make_carbonate_thickness_video.py --max-age 80 --framerate 10

Requires: gplately, pygplates, cartopy, xarray, matplotlib, and ffmpeg on PATH.
Outputs (in <figures>/videos/):
    carbonate_thickness_back_in_time.mp4      (present -> past)
    carbonate_thickness_forward_in_time.mp4   (past -> present)
"""
import argparse, os, sys
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _ccd_video_common as V

# ---- path resolution -------------------------------------------------------
# Works both in the working tree (this script in Paper/, the workflow in a
# sibling CCD_workflow_clean/) and in the public repo (this script in
# figure_scripts/, the workflow at the repo root). CW = workflow root,
# FIGDIR = where figures and the small figure-input curves live.
_HERE = Path(__file__).resolve().parent
def _workflow_root():
    for _p in (_HERE.parent / "CCD_workflow_clean", _HERE.parent, _HERE):
        if (_p / "steps").is_dir() and (_p / "ccdworkflow").is_dir():
            return _p
    raise SystemExit("cannot locate the CCD workflow root (expected steps/ + ccdworkflow/)")
CW = _workflow_root()
FIGDIR = _HERE / "Figures" if (_HERE / "Figures").is_dir() else CW / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR = CW / "steps/step9_carbonate_volume_analysis/input/Alfonso_etal_2024_modClennettMuller"
CPT       = CW / "data/carbonate_thickness_roma_pale.cpt"
GRID_DIR  = CW / "steps/step8_carbonate_sediment_thickness/carbonate_sed_thickness_DM2026"
GRID_FMT  = str(GRID_DIR / "compacted_sediment_thickness_0.25_{t}.nc")
OUT       = FIGDIR / "videos"; OUT.mkdir(parents=True, exist_ok=True)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-age", type=int, default=170)
    ap.add_argument("--time-step", type=int, default=1)
    ap.add_argument("--framerate", type=int, default=8)
    ap.add_argument("--force", action="store_true", help="re-render existing frames")
    a = ap.parse_args()
    cmap, norm = V.gmt_cpt(str(CPT))
    times = list(range(0, a.max_age + 1, a.time_step))
    print(f"carbonate-thickness video: {len(times)} frames (0-{a.max_age} Ma / {a.time_step} Myr)")
    V.render_video(times, GRID_FMT, str(OUT), "carbonate_thickness",
                   "Compacted carbonate sediment thickness (m)", str(MODEL_DIR),
                   cmap, norm=norm,
                   # Every class the same width, and every boundary labelled, as in
                   # Figure 3: the ten classes below 50 m are the ones to tell apart.
                   cbar_ticks=[0,5,10,15,20,25,30,35,40,45,50,100,160,210,270],
                   cbar_extend="max", cbar_spacing="uniform",
                   framerate=a.framerate, force=a.force)

if __name__ == "__main__":
    main()
