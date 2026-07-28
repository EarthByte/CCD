# Changelog

Divergence of this local copy from the upstream `EarthByte/CarbonateSedimentThickness`
repository (Dutkiewicz et al., 2019, *Geology*, doi:10.1130/G45424.1). Each
entry lists what changed, why, and where to look in the code.

The intent is to merge these changes upstream — see `email_to_developer.md`
for the cover note to John Cannon.

---

## [Unreleased] — local additions for the Dutkiewicz & Müller (in prep) CCD paper

### Added — latitude- and time-dependent correction (`latitude_factor.py`)

New standalone module providing `latitude_factor(lat_deg, time_ma=0.0, ...)`,
a multiplier in [0, 1] applied to the per-step carbonate sedimentation rate
inside `calc_carbonate_decompacted_sediment_thickness`. The factor is the
product of a spatial taper and a temporal blend factor:

- **Spatial half-cosine taper in absolute palaeolatitude.** Returns 1.0
  equatorward of `HIGH_LATITUDE_FULL_PRESERVATION_DEG` (default 60°),
  smoothly decreases to 0.0 at `HIGH_LATITUDE_ZERO_ACCUMULATION_DEG`
  (default 80°). Calibrated against modern sediment-trap CaCO₃ flux
  compilations (Honjo et al. 2008), the modern CaCO₃ atlas
  (Archer 1996), and the absence of pelagic carbonate in the modern
  deep Arctic (Stein 2008; Jutterström & Anderson 2016).

- **Temporal half-cosine ramp in palaeotime.** Disabled (factor 1.0
  everywhere) at and before `GLACIATION_ONSET_MA` (default 34 Ma),
  fully active at and after `FULL_CORRECTION_MA` (default 23 Ma),
  smoothly ramped in between. The 34 Ma anchor is the EOT-1 / Oi-1
  Antarctic glaciation onset (Coxall et al. 2005; Miller et al. 2008);
  the 23 Ma anchor is the end of the Oligocene pCO₂ drawdown
  reconstructed by the CenCO2PIP consortium (Hönisch et al. 2023).
  Direct evidence that pelagic carbonate was being deposited at >55° S
  palaeolatitude through the pre-34 Ma greenhouse comes from Maud Rise
  ODP 689/690 (Thomas 1990), the Mentelle Basin IODP Expedition 369
  sites U1513–U1516 (Huber et al. 2019), and the Campbell Plateau /
  Otway Basin records that localise the carbonate shutdown to the EOT
  (Gallagher et al. 2020). Without the temporal ramp a
  modern-calibrated restriction would over-strip the geological record
  in those greenhouse intervals.

Both tapers are C¹ continuous so the time-derivative of the correction
is zero at the four threshold ages/latitudes — no kinks propagate into
time-derivative diagnostics (e.g. plate-influx rates).

Full DOIs for all cited references are in the `latitude_factor.py`
module docstring.

### Changed — `carbonate_sediment_thickness.py`

- Added `from latitude_factor import latitude_factor` to the imports.
- Inside `calc_carbonate_decompacted_sediment_thickness` (around the
  per-step rate calculation, originally lines ~519–525), the
  `carbonate_decompacted_sediment_rate` is now multiplied by
  `latitude_factor(paleo_lat, reconstruction_time)` where `paleo_lat`
  comes from `reconstructed_points[index].to_lat_lon()`. This is the
  *paleo*-latitude of the reconstructed point at the current
  reconstruction time, so crust which originated at low latitude still
  receives credit for carbonate accumulated before it drifted poleward.

### Changed — `carbonate_sediment_thickness.py` (gridline registration made explicit)

In `write_grid_file_from_xyz` the `gmt nearneighbor` call now carries an
explicit `-rg` flag (gridline registration). The default is already
gridline, but stating it explicitly:

1. Documents the intent in the source rather than relying on a default
   that could change in a future GMT release; and
2. Makes the registration unambiguous to downstream tools, which
   eliminates the `Longitude range too small; geographic boundary
   condition changed to natural` warning that `grdimage` emits when a
   global grid's data span is not 360°.

Downstream consumers (the video and panel scripts in
`8_carbonate_thickness_analysis/`) assume this gridline convention.

### Sensitivity-test parameters

Both the spatial cut-off and the temporal ramp are exposed as
module-level constants in `latitude_factor.py` so they can be varied
without editing the inner loop:

| Parameter                              | Default | Defensible sweep        |
| -------------------------------------- | ------- | ----------------------- |
| `HIGH_LATITUDE_FULL_PRESERVATION_DEG`  | 60°     | 50° / **60°** / 65°      |
| `HIGH_LATITUDE_ZERO_ACCUMULATION_DEG`  | 80°     | 70° / **80°** / 85°      |
| `GLACIATION_ONSET_MA`                  | 34.0    | held fixed at EOT-1/Oi-1 |
| `FULL_CORRECTION_MA`                   | 23.0    | 33 / **23** / 16         |

### Added — `CHANGELOG.md`

Repo-level documentation tracking the divergence from the upstream workflow.

---

## Upstream baseline

All other files are unchanged from upstream
[EarthByte/CarbonateSedimentThickness](https://github.com/EarthByte/CarbonateSedimentThickness)
as of the 2019_v2 topology model release.
