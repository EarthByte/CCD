# CCD — global carbonate compensation depth evolution over the last 170 million years and its drivers

Code, derived results and figures for:

> Dutkiewicz, A., and Müller, R.D., *Global carbonate compensation depth evolution
> over the last 170 million years and its drivers* (submitted). EarthByte Group,
> School of Geosciences, The University of Sydney.

The repository reconstructs an area-weighted global carbonate compensation depth
(CCD) from Atlantic, Pacific and Indian Ocean records, calibrates it against
long-term sea level, extends it to 170 Ma, converts it into global carbonate
sediment thickness, and drives a solid-Earth carbon-degassing model with the
result. It contains everything needed to reproduce the analysis except the large
input and output grids, which are archived on Zenodo (see
[Data availability](#data-availability)).

---

## Quick start

```bash
git clone https://github.com/EarthByte/CCD.git
cd CCD
pip install -r requirements.txt

python run_ccd_core.py          # steps 1–7: timescale normalisation → global CCD → 170 Ma hybrid
```

`run_ccd_core.py` is self-contained and deterministic: every input it needs is in
`data/`, and it regenerates every step-1–7 output and diagnostic figure. It takes
a few minutes on a laptop.

The downstream stages need a local geoscience stack and the Zenodo grids:

```bash
./run_carbon_cycle.sh           # steps 8–10: thickness grids → volumes → CO2 degassing
```

`run_carbon_cycle.sh` re-runs `run_ccd_core.py` first, so the carbon stages always
start from the current CCD and nothing stale propagates. See
`pipeline_carbon/README.md` for its options (`--skip-ccd`, `--only-step9`,
`CO2_START_FROM=…` to resume a notebook sequence).

---

## Workflow

| Step | Folder | What it does |
|-----:|--------|--------------|
| 1 | `steps/step1_timescale_conversion` | Normalise every CCD and sea-level curve to GTS2020, then resample onto a common grid |
| 2 | `steps/step2_ocean_basin_areas` | Reconstructed ocean-basin area fractions through time (notebook; result stored in `data/ccds_and_oceanbasin_fractions.xlsx`) |
| 3 | `steps/step3_global_ccd_synthesis` | Area-weighted global CCD synthesis and inter-basin dispersion, 0–52 Ma |
| 4 | `steps/step4_sealevel_envelope` | Hybrid Miller et al. (2024) / Haq et al. (1987) sea-level peak-following envelope |
| 5 | `steps/step5_ccd_lowpass_filter` | Butterworth low-pass of the global CCD to match the sea-level bandwidth |
| 6 | `steps/step6_sealevel_ccd_regression` | Reduced-major-axis sea-level–CCD calibration, prediction and diagnostics |
| 7 | `steps/step7_planktogenic_hybrid_ccd` | Planktogenic correction and splice → hybrid CCD, 0–170 Ma |
| 8 | `steps/step8_carbonate_sediment_thickness` | Compacted carbonate sediment thickness grids from the new CCD (min / mean / max sedimentation-rate scenarios) |
| 9 | `steps/step9_carbonate_volume_analysis` | Global carbonate volume and area statistics; DSDP/ODP ground-truth at 38 basement sites |
| 10 | `steps/step10_carbon_cycle_degassing` | Deep-Earth carbon-cycle notebooks: subducted carbon → degassing → atmospheric flux |

Steps 1 and 3–7 are pure Python (numpy / scipy / pandas / matplotlib). Step 2 is a
notebook. Steps 8–10 additionally need `pygplates`, `gplately`, GMT/`pygmt` and
`jupyter`, plus the Zenodo grids.

### Timescales

The published source curves sit on three different geomagnetic polarity
timescales: GTS2020 (Ogg, 2020) for the Miller et al. (2024) sea level and the
Dalvand et al. (2025) Pacific and Indian CCDs, GTS2012 (Ogg, 2012) for the
Atlantic CCD and the Haq sea-level compilations, and Cande & Kent (1995) for the
Pälike et al. (2012) record that extends the Pacific CCD from 35.5 to 52 Ma. **Step 1 is
the single place where these age models are reconciled.** Everything downstream
reads GTS2020 curves from `steps/step1_timescale_conversion/outputs/`. Only curves
that enter the analysis are converted; published curves shown purely for comparison,
including the Haq (1987) pair used to validate the sea-level envelope, keep their own
age models. Segments
already on GTS2020 pass through untouched, splices are re-imposed after conversion,
and the converted series are resampled onto a 1 Myr grid for the CCDs and 0.1 Myr
for sea level. See `steps/step1_timescale_conversion/README.md`.

Each step writes to its own `outputs/` (or `output/`) folder. Those outputs are
committed here so the results can be inspected without rerunning the workflow.

**Notebooks are committed with their outputs cleared** to keep the repository
small. Their executed results are included as CSV under
`steps/step10_carbon_cycle_degassing/Alfonso_etal_2024_DM26/Outputs/Notebook*/csv/`.

---

## Figures

`figure_scripts/` holds the scripts that build the published figures and the two
supplementary animations; `figures/` holds the rendered output (PDF + PNG), the
hybrid CCD curve and the carbonate-budget table that panel 3D plots.

Figures 1, 2 and 4 are drawn with pyGMT at the journal's column widths, so the
type sizes on the page are the sizes the reader sees; Figure 3 and the
supplementary figures are matplotlib. All of them are set in Helvetica.

| Figure | Script | Output |
|--------|--------|--------|
| Fig. 1 — regional and area-weighted global CCD, 0–52 Ma | `figure_scripts/make_fig1_gmt.py` | `figures/Fig1_regional_global_ccd_gmt.*` |
| Fig. 2 — long-term sea level with the CCD calibration inset, and the hybrid CCD to 170 Ma | `figure_scripts/make_fig2_gmt.py` | `figures/Fig2_combined_gmt.*` |
| Fig. 3 — reconstructed carbonate thickness at 0, 34 and 115 Ma, with the carbonate budget | `figure_scripts/make_fig3_maps.py` (panel D via `carbonate_budget.py`) | `figures/Fig3_carbonate_thickness_maps.*`, `figures/Fig3_carbonate_budget.csv` |
| Fig. 4 — the CCD against the solid-Earth carbon fluxes, their correlations, and the variance each explains | `figure_scripts/make_fig4_gmt.py` | `figures/Fig4_combined_gmt.*` |
| Fig. S1 — regression-method sensitivity | `steps/step6_sealevel_ccd_regression/regression_diagnostics.py` | `figures/FigS1_regression_method_sensitivity.*` |
| Fig. S2 — modelled vs observed present-day thickness | `steps/step9_carbonate_volume_analysis/08_ground_truth_carbonate_thickness.py` | `figures/FigS2_modelled_vs_observed_thickness.*` |
| Fig. S3 — residual distribution at 38 sites | `steps/step9_carbonate_volume_analysis/08_ground_truth_carbonate_thickness.py` | `figures/FigS3_thickness_residual_histogram.*` |
| Fig. S4 — residuals by ocean basin | `steps/step9_carbonate_volume_analysis/08_ground_truth_carbonate_thickness.py` | `figures/FigS4_thickness_residuals_by_basin.*` |
| Fig. S5 — net-outflux end-member scenarios | `steps/step10_carbon_cycle_degassing/standalone_plots` | `figures/FigS5_net_outflux_endmembers.*` |

`figures/Fig1_regional_global_ccd.*` is the same panel rendered by
`steps/step3_global_ccd_synthesis/global_ccd_synthesis.py` as that step's own
diagnostic.

`figure_scripts/attribution_closure_analysis.py` performs the attribution
(correlations, AR1-adjusted significance, incremental R², first differences) and
writes `figure_scripts/attribution_matched_series.csv` and
`figure_scripts/component_correlations.csv`. Figure 4 reads both, so the panels and
the numbers quoted in the text cannot drift apart.

The figure scripts carry their own layout checks — text-collision detection, minimum
type size against the journal's floor, line weights and the typefaces embedded in the
PDF — and report them when they run.

### Animations

| Movie | Script |
|-------|--------|
| Movie S1 — global paleobathymetry, 170–0 Ma | `figure_scripts/make_paleobathymetry_video.py` |
| Movie S2 — compacted carbonate sediment thickness, 170–0 Ma | `figure_scripts/make_carbonate_thickness_video.py` |

Both share `figure_scripts/_ccd_video_common.py` and render 171 frames plus MP4s
into `figures/videos/` (git-ignored). The rendered movies are on Zenodo. Options:

```bash
python figure_scripts/make_carbonate_thickness_video.py --time-step 5   # quick preview
python figure_scripts/make_paleobathymetry_video.py --force             # re-render all frames
```

The figure and video scripts locate the workflow root and the figure folder
automatically, so they run unchanged from this repository or from the authors'
working tree.

---

## Supplementary datasets

The three datasets cited by the supplement are in `supplementary/`:

| Dataset | File | Contents |
|---------|------|----------|
| S1 | `supplementary/Supplementary_Dataset_S1_CCD.xlsx` | the global and regional CCD curves and the reconstructed ocean-basin area fractions |
| S2 | `supplementary/Supplementary_Dataset_S2_carbon_model.xlsx` | the deep-Earth carbon-cycle model results: plate influx, atmospheric influx by source, and net outflux |
| S3 | `supplementary/Supplementary_Dataset_S3_carbonate_budget.csv` | the carbonate budget plotted in Figure 3D — net volume gain and seafloor area above the CCD |

S1 and S2 are rebuilt from the current workflow outputs by
`figure_scripts/make_supplementary_datasets.py`; every sheet records the file its
numbers came from, and each workbook carries a provenance sheet that marks any
source predating the current CCD. S3 is written by `carbonate_budget.py` at the
same time as the figure panel, from the same arrays.

---

## Data availability

Small inputs — the regional CCD curves, the sea-level records, the reference CCDs
of Boss & Wilkinson (1991) and Delaney & Boyle (1988), the ocean-basin area
fractions and the present-day carbonate-thickness validation grids — are in
`data/` and are all that `run_ccd_core.py` needs.

The large grids are archived on Zenodo:

> **Zenodo archive:** *DOI to be assigned* — `https://doi.org/10.5281/zenodo.XXXXXXX`

It contains:

- **Paleobathymetry grids**, 170–0 Ma at 1 Myr, on the Alfonso et al. (2025) plate model
- **Compacted carbonate sediment thickness grids**, 170–0 Ma at 1 Myr (minimum, mean and maximum sedimentation-rate scenarios)
- **Decompacted carbonate sediment thickness grids**, same coverage
- the plate model, continental masks and reservoir grids used by steps 8–10
- the rendered supplementary animations

Download the archive and place its contents at the paths given in
`pipeline_carbon/config.sh` before running steps 8–10.

---

## Citation

If you use this code or the derived CCD record, please cite the paper and this
repository (see `CITATION.cff`), and the Zenodo archive for the grids.

## Licence

Code is released under the MIT Licence (`LICENSE`). Figures, derived data tables
and the reconstructed CCD record are released under CC BY 4.0
(`LICENSE-FIGURES`). Third-party inputs retain their original licences; the
carbonate-sediment-thickness workflow in
`steps/step8_carbonate_sediment_thickness/` carries its own upstream `LICENSE`.
