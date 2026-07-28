#!/usr/bin/env bash
# ============================================================================
# run_all.sh -- run the full carbonate-thickness comparison pipeline.
#
# Pipeline (dependency order):
#   02  per-time difference stats + 4 picked times
#   03  3-panel boxplot figure
#   04  regional (basin x lat-band) heatmap + per-picked-time tables
#   06  2x2 panel of delta maps at the picked times
#   08  ground-truth of present-day carbonate thickness vs DSDP/ODP well data
#   05  MS Word writeup
#
# Videos (01) are opt-in -- pass --videos to render them too.
#
# Usage:
#   ./run_all.sh                 # stats + figs + writeup (no videos)
#   ./run_all.sh --videos        # everything, including the three MP4s
#   ./run_all.sh --quick-videos  # videos at 5 Myr cadence (smoke test)
# ============================================================================
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

DO_VIDEOS=0
QUICK_VIDEOS=0
for a in "$@"; do
  case "$a" in
    --videos)       DO_VIDEOS=1 ;;
    --quick-videos) DO_VIDEOS=1; QUICK_VIDEOS=1 ;;
    *) echo "Unknown arg: $a" >&2; exit 2 ;;
  esac
done

echo "==> 02: per-time difference statistics + picked times"
python3 02_difference_stats.py

echo
echo "==> 03: 5 Myr-binned boxplots"
python3 03_boxplots.py

echo
echo "==> 04: regional (basin x lat-band) characterisation"
python3 04_region_characterisation.py

echo
echo "==> 06: 2x2 panel of delta maps at the picked times"
python3 06_delta_panel.py

echo
echo "==> 08: ground-truth vs DSDP/ODP well data (map, histogram, boxplots)"
python3 08_ground_truth_carbonate_thickness.py

echo
echo "==> 05: MS Word writeup"
python3 05_build_writeup.py

if [ "$DO_VIDEOS" = "1" ]; then
  echo
  if [ "$QUICK_VIDEOS" = "1" ]; then
    echo "==> 01: rendering videos at 5 Myr cadence (smoke test)"
    python3 01_render_videos.py --cadence 5 --force
  else
    echo "==> 01: rendering all three videos at 1 Myr cadence"
    python3 01_render_videos.py
  fi
else
  echo
  echo "(skipping videos -- pass --videos or --quick-videos to render them)"
fi

echo
echo "Done.  Outputs in $(pwd)/output/"
