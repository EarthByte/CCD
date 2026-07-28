#!/usr/bin/env python3
"""
Supplementary Movie S? -- global paleobathymetry, 0-170 Ma.

Mollweide projection, coloured with the `nuuk` scientific colour map over
-5000..0 m (deeper clipped to the deepest colour, subaerial clipped to the
shallowest), matching the pyBacktrack paleobathymetry-video template. Same
GPlately overlays as the carbonate video: reconstructed plate boundaries
(black + white halo), subduction teeth, gray continental polygons + thin
coastlines. NaN cells are gray (matching the template's NAN colour). The
paleobathymetry grids are already in Alfonso et al. (2025) palaeo-coordinates.

Usage:
    python3 make_paleobathymetry_video.py                 # 0-170 Ma, 1 Myr, 8 fps
    python3 make_paleobathymetry_video.py --time-step 5   # quick preview

Requires: gplately, pygplates, cartopy, xarray, matplotlib, cmcrameri, ffmpeg.
Outputs (in <figures>/videos/):
    paleobathymetry_back_in_time.mp4      (present -> past)
    paleobathymetry_forward_in_time.mp4   (past -> present)
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
GRID_DIR  = CW / "steps/step8_carbonate_sediment_thickness/input_grids/Alfonso2024_pybacktrack_merged_paleobathymetry"
GRID_FMT  = str(GRID_DIR / "paleobathymetry_{t}Ma.nc")
OUT       = FIGDIR / "videos"; OUT.mkdir(parents=True, exist_ok=True)

def _bathy_cmap():
    try:
        from cmcrameri import cm as cmc
        cmap = cmc.nuuk.copy()
    except Exception:
        import matplotlib.pyplot as plt
        print("  [note] cmcrameri not found; falling back to 'YlGnBu_r'")
        cmap = plt.get_cmap("YlGnBu_r").copy()
    cmap.set_bad("0.784")     # 200/200/200, matches template NAN_COLOR
    return cmap

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-age", type=int, default=170)
    ap.add_argument("--time-step", type=int, default=1)
    ap.add_argument("--framerate", type=int, default=8)
    ap.add_argument("--force", action="store_true", help="re-render existing frames")
    a = ap.parse_args()
    cmap = _bathy_cmap()
    times = list(range(0, a.max_age + 1, a.time_step))
    print(f"paleobathymetry video: {len(times)} frames (0-{a.max_age} Ma / {a.time_step} Myr)")
    V.render_video(times, GRID_FMT, str(OUT), "paleobathymetry",
                   "Paleobathymetry (m)", str(MODEL_DIR),
                   cmap, vmin=-5000, vmax=0,
                   cbar_ticks=[-5000,-4000,-3000,-2000,-1000,0],
                   # only the deep end is open-ended: there are no land elevations
                   # in these grids, so no arrow at the 0 m end of the bar
                   cbar_extend="min", framerate=a.framerate, force=a.force)

if __name__ == "__main__":
    main()
