# Step 1 — timescale normalisation to GTS2020

The CCD and sea-level records this study combines were published on three
different geomagnetic polarity timescales. Splicing or regressing them without a
common age model mixes those timescales into the result, so every curve is mapped
onto GTS2020 (Ogg, 2020) here, before anything else runs.

Run it with

```
python steps/step1_timescale_conversion/convert_to_gts2020.py
```

or simply `python run_ccd_core.py`, which runs it first.

## What is converted

| Curve | Source | Original timescale | Action |
|-------|--------|--------------------|--------|
| Atlantic CCD | Dutkiewicz & Müller (2022) | GTS2012 | converted |
| Pacific CCD, 0–35.5 Ma | Dalvand et al. (2025) | GTS2020 | passed through |
| Pacific CCD, 35.5–52 Ma | Pälike et al. (2012) | Cande & Kent (1995) | converted |
| Indian CCD | Dalvand et al. (2025) | GTS2020 | passed through |
| Hybrid short-term sea level, 0–66.6 Ma | Miller et al. (2024) | GTS2020 | passed through |
| Hybrid short-term sea level, >66.6 Ma | Haq et al. | GTS2012 | converted |

A segment that is already on GTS2020 is passed through untouched, so no spurious
interpolation is introduced where none is needed.

## What is deliberately not converted

Only curves that enter the analysis are converted. Curves used solely for comparison
keep their published age models:

- the standalone Haq (1987) short-term and long-term sea-level curves, which are the
  step-4 validation pair. The envelope operator is applied to one and checked against
  the other, and both are on GTS2012 as published, so the check is already
  self-consistent. Remapping them changes the reported statistic through resampling
  alone: r 0.900 to 0.898, RMSE 35.5 to 36.4 m, MAE 27.1 to 27.6 m, with no change to
  anything the comparison is testing.
- the Boss & Wilkinson (1991) and Delaney & Boyle reference CCDs plotted in Fig. 2b.
  Neither of their timescales is carried by `timescales.txt`, so neither could be
  converted in any case.

Note that the Haq compilation still reaches the analysis, as the pre-66.6 Ma half of
the hybrid short-term sea-level record, and that half is converted. It is the
standalone published curves that are left alone.

## How the ages are mapped

Ages are remapped by piecewise-linear interpolation between the polarity-chron
boundaries that the two timescales share, using the EarthByte
`timescale_conversion.py` tool and its `timescales.txt` chron table. Both files
are copied into this folder so the workflow is self-contained; the originals live
in `scripts/1_timescale_conversion/timescale_conversion/`. The table carries
Gee & Kent (2007), Cande & Kent (1995), Gradstein (2005), GTS2012 and GTS2020.

Two details matter in practice.

**Splices are re-imposed after conversion.** A splice age between two records is
defined in the original frame. Converting the older segment can move its youngest
samples back across that boundary, where they would interleave with the younger
record and put a sawtooth into the spliced curve. Any converted sample that lands
at or below the highest age already accepted is therefore dropped, and the regrid
bridges the short gap linearly. For the Miller/Haq sea-level splice at 66.6 Ma
this drops a single sample and leaves a 0.08 Myr gap.

**The remap leaves the samples unevenly spaced**, so each series is put back onto
a uniform grid: 1 Myr for the CCD curves, 0.1 Myr for the sea-level series. The
sea-level grid has to stay dense because step 4 runs a rolling-quantile envelope
over it, which a 1 Myr grid would destroy.

## Inputs

```
data/ccds_and_oceanbasin_fractions.xlsx            regional CCD columns 0-5
data/sealevel/Miller_Haq_SeaLevel_ShortTerm_hybrid.tsv
```

The regional CCDs are read from the spreadsheet rather than from the archived
0.5 Myr files, because the spreadsheet columns are the curves that produced the
published synthesis: for the Pacific and Indian the two agree except in the
youngest few Myr, but the Atlantic column is a later revision differing by up to
741 m.

## Outputs (`outputs/`)

```
ATL_CCD_GTS2020_1my.txt                 Age_Ma, CCD_m       0-65 Ma
PAC_CCD_GTS2020_1my.txt                                     0-52 Ma
IND_CCD_GTS2020_1my.txt                                     0-23 Ma
sealevel_shortterm_hybrid_GTS2020.txt   Age_Ma, Sea_level_m 0-205 Ma, 0.1 Myr
timescale_conversion_report.txt         what moved, and by how much
timescales.log                          chron table as parsed, from the tool
```

## Size of the correction

Over the 0–52 Ma calibration interval the conversion moves the Atlantic CCD by up
to 139 m (rms 39 m) and the Pacific by up to 153 m (rms 36 m); the Indian is
unchanged, being GTS2020 already. In the Cretaceous the sea-level ages move by up
to 4.5 Myr, and by 1.3 Myr on average. `timescale_conversion_report.txt` records
the shift for every segment.

## Who reads these outputs

```
step3_global_ccd_synthesis   the three regional CCD curves
step4_sealevel_envelope      the hybrid short-term sea-level record
Paper/make_fig2_combined.py  the hybrid short-term sea-level record
```

