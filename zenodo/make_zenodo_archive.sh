#!/usr/bin/env bash
# ============================================================================
# make_zenodo_archive.sh — assemble the Zenodo grid archive for the CCD paper.
#
#   ./zenodo/make_zenodo_archive.sh [output_dir]
#   ./zenodo/make_zenodo_archive.sh --check       # verify the sources, archive nothing
#
# Default output_dir is ../CCD_zenodo_archive, i.e. beside the repository and
# outside it, so a multi-gigabyte archive is never a candidate for committing.
#
# Each component is tarred straight from its source folder - nothing is copied
# first - then hashed. The script counts the files it is about to archive and
# stops if a component is missing or short, so a half-finished run cannot be
# uploaded as if it were complete.
# ============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
CHECK_ONLY=0
if [ "${1:-}" = "--check" ]; then CHECK_ONLY=1; shift; fi
OUT="${1:-$(cd "$REPO/.." && pwd)/CCD_zenodo_archive}"

# Where the grids live. The heavy grids are git-ignored, so in the authors' working
# tree they sit in the sibling workflow folder rather than inside this repository.
# Look in the repository first, then in the sibling, then give up and say so.
# Setting STEPS_ROOT overrides the search, and matches pipeline_carbon/config.sh.
if [ -z "${STEPS_ROOT:-}" ]; then
  for cand in "$REPO/steps" "$(cd "$REPO/.." && pwd)/CCD_workflow_clean/steps"; do
    if [ -d "$cand/step8_carbonate_sediment_thickness/carbonate_sed_thickness_DM2026" ]; then
      STEPS_ROOT="$cand"; break
    fi
  done
fi
STEPS_ROOT="${STEPS_ROOT:-$REPO/steps}"
echo "Grids: $STEPS_ROOT"
S8="$STEPS_ROOT/step8_carbonate_sediment_thickness"
S9="$STEPS_ROOT/step9_carbonate_volume_analysis"
S10="$STEPS_ROOT/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26"
# The rendered movies likewise live wherever the figure scripts wrote them.
if [ -z "${VIDEOS:-}" ]; then
  for cand in "$REPO/figures/videos" "$(cd "$REPO/.." && pwd)/Paper/Figures/videos" \
              "$(cd "$REPO/.." && pwd)/CCD_workflow_clean/figures/videos"; do
    if [ -f "$cand/carbonate_thickness_back_in_time.mp4" ]; then VIDEOS="$cand"; break; fi
  done
fi
VIDEOS="${VIDEOS:-$REPO/figures/videos}"
echo "Movies: $VIDEOS"

# name | parent directory | folder to archive | expected file count (0 = don't check)
COMPONENTS=(
  "carbonate_sediment_thickness_mean|$S8|carbonate_sed_thickness_DM2026|513"
  "carbonate_sediment_thickness_min|$S8|carbonate_sed_thickness_min_DM2026|513"
  "carbonate_sediment_thickness_max|$S8|carbonate_sed_thickness_max_DM2026|513"
  "paleobathymetry|$S8/input_grids|Alfonso2024_pybacktrack_merged_paleobathymetry|171"
  "continental_masks|$S10/Grids/InputGrids|ContinentalMasks|171"
  "plate_model|$S9/input|Alfonso_etal_2024_modClennettMuller|0"
  "present_day_validation_grids|$REPO/data|carbonate_thickness|0"
)

if [ "$CHECK_ONLY" = "0" ]; then # repack_grids.py writes scaled-integer copies into <step8>/repacked/. When they are
# there they are what gets archived: same values to within 0.1 m, less than half the
# size. Delete that folder to fall back to the float32 originals.
REPACKED="$S8/repacked"
resolve_parent() {   # $1 = default parent, $2 = folder name
  if [ -d "$REPACKED/$2" ]; then printf '%s' "$REPACKED"; else printf '%s' "$1"; fi
}

mkdir -p "$OUT"; echo "Writing the archive to $OUT"; fi
echo

missing=()
for spec in "${COMPONENTS[@]}"; do
  IFS='|' read -r name parent folder expected <<< "$spec"
  parent="$(resolve_parent "$parent" "$folder")"
  src="$parent/$folder"
  if [ ! -d "$src" ]; then missing+=("$name: no such folder $src"); continue; fi
  n=$(find "$src" -type f ! -name '.DS_Store' | wc -l | tr -d ' ')
  if [ "$expected" -gt 0 ] && [ "$n" -ne "$expected" ]; then
    missing+=("$name: $n files, expected $expected")
  fi
done

# The animations are four named files rather than a folder.
MOVIES=(carbonate_thickness_back_in_time.mp4 carbonate_thickness_forward_in_time.mp4
        paleobathymetry_back_in_time.mp4 paleobathymetry_forward_in_time.mp4)
for m in "${MOVIES[@]}"; do
  [ -f "$VIDEOS/$m" ] || missing+=("animations: $VIDEOS/$m not found")
done

if [ ${#missing[@]} -gt 0 ]; then
  echo "Nothing archived. Fix these first:" >&2
  printf '  %s\n' "${missing[@]}" >&2
  exit 1
fi
echo "All components present."
if [ "$CHECK_ONLY" = "1" ]; then
  for spec in "${COMPONENTS[@]}"; do
    IFS='|' read -r name parent folder expected <<< "$spec"
    parent="$(resolve_parent "$parent" "$folder")"
    printf '  %-32s %6s  %s\n' "$name" "$(du -sh "$parent/$folder" | cut -f1)" "$parent/$folder"
  done
  printf '  %-32s %6s  %s\n' animations "$(du -ch "${MOVIES[@]/#/$VIDEOS/}" | tail -1 | cut -f1)" "$VIDEOS"
  exit 0
fi

for spec in "${COMPONENTS[@]}"; do
  IFS='|' read -r name parent folder expected <<< "$spec"
  parent="$(resolve_parent "$parent" "$folder")"
  echo "  $name  <-  $parent/$folder"
  tar --exclude='.DS_Store' -czf "$OUT/$name.tar.gz" -C "$parent" "$folder"
done

echo "  animations  <-  $VIDEOS"
tar -czf "$OUT/animations.tar.gz" -C "$VIDEOS" "${MOVIES[@]}"

cp "$HERE/README_zenodo.md" "$OUT/README.md"

echo
echo "Hashing..."
(
  cd "$OUT"
  : > MANIFEST.txt
  for f in *.tar.gz README.md; do
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$f" >> MANIFEST.txt
    else shasum -a 256 "$f" >> MANIFEST.txt; fi
  done
)

echo
echo "Archive contents:"
ls -lh "$OUT" | tail -n +2
echo
echo "Total: $(du -sh "$OUT" | cut -f1)"
echo
echo "Upload every file in $OUT to the Zenodo deposition, README.md included."
