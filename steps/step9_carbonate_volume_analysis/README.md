# 8_carbonate_thickness_analysis

Comparison of global compacted-carbonate-sediment-thickness grids produced
by two alternative CCD curves (Dutkiewicz & Müller 2026, in prep, and Boss
& Wilkinson 1991), both built on the same Alfonso et al. (2024)
pyBacktrack paleobathymetry.

## Pipeline

```
config.py                  shared paths, time list, region helpers
paleo_coastlines.py        Alfonso 2024 coastline reconstruction helper
01_render_videos.py        3 MP4 videos (DM2026, BW1991, DM2026 − BW1991)
                           Winkel-Tripel projection, reconstructed coastlines
02_difference_stats.py     per-time stats CSV + 4 picked times sidecar
03_boxplots.py             2-panel fig: (a) DM/BW mean ± σ time-series
                           (b) signed delta boxplot, 5 Myr bins
04_region_characterisation.py
                           regional (basin × lat-band) heatmap + tables
06_delta_panel.py          2 × 2 global delta maps at the 4 picked times
05_build_writeup.py        assembles output/carbonate_thickness_comparison.docx
```

Dependency order: 02 must run before 03, 04, 06 and 05; 03, 04 and 06 must
run before 05. `01_render_videos.py` is independent.

## Quick start

```bash
cd CCD_workflow_clean/steps_carbon/step8_analysis

# Stats + picked times first (under a minute)
python 02_difference_stats.py

# Figures (a few minutes each)
python 03_boxplots.py
python 04_region_characterisation.py
python 06_delta_panel.py

# Word writeup (a few seconds)
python 05_build_writeup.py

# Videos (~8 min each, opt-in)
python 01_render_videos.py                     # all three
python 01_render_videos.py --modes dm2026      # just one
python 01_render_videos.py --cadence 5 --modes delta --force   # quick smoke test
```

## Outputs

```
output/
├── carbonate_thickness_comparison.docx     <- the writeup
├── carbonate_dm2026_150-0Ma.mp4            <- DM2026 thickness video
├── carbonate_bw1991_150-0Ma.mp4            <- BW1991 thickness video
├── carbonate_delta_150-0Ma.mp4             <- DM2026 - BW1991 video
├── figures/
│   ├── 03_thickness_and_delta.{png,pdf}    <- thickness mean ± σ + delta boxplot
│   ├── 04_region_heatmap.{png,pdf}
│   └── 06_delta_panel.{png,pdf}            <- 2 × 2 delta maps at picked times
└── stats/
    ├── difference_stats.csv                <- per-1Myr delta stats
    ├── picked_times.txt                    <- 4 NMS-picked ages
    ├── thickness_mean_std.csv              <- per-1Myr DM/BW mean & std
    ├── boxplot_delta_summary.csv           <- per-5Myr-bin delta stats
    ├── regional_mean_abs.csv               <- wide (time × basin × lat-band) matrix
    └── picked_time_region_tables.csv       <- top-3 regions per picked time
```

## Key conventions

- **Difference sign**: positive Δ = DM2026 produces more carbonate than BW1991.
- **Time range**: 0–150 Ma (the grids exist to 170 Ma but pre-150 Ma carbonate
  area is small).
- **Common mask**: all summary statistics use only cells where BOTH products
  produce carbonate at that time. Cells where one product is NaN (above CCD,
  on land, etc.) and the other isn't are excluded from the comparison.
- **Picked times**: greedy non-maximum-suppression on per-time mean |Δ|,
  exclusion radius 25 Myr, n = 4.
- **Basins**: present-day Pacific / Atlantic / Indian / Arctic / Southern via
  simple lon/lat cuts in `config.basin_of`. No paleo-reconstruction of the
  cells.

## Software prerequisites

```bash
# pyGMT for the videos
conda install -c conda-forge pygmt
# ffmpeg for stitching MP4s
conda install -c conda-forge ffmpeg
# python-docx for the writeup
pip install python-docx
# xarray + netCDF4 for the grid I/O (probably already present)
conda install -c conda-forge xarray netcdf4
```
