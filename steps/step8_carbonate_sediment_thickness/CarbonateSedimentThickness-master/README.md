# Calculate carbonate sediment thickness

Generate carbonate sediment thickness grids from age and bathymetry grids over the time range 0-170Ma (in 1My increments). Pelagic carbonate sediments did not exist before 170Ma.

## Dependencies

- [GMT](https://www.generic-mapping-tools.org/download/) (and make sure the 'gmt' executable is in the PATH).
- [PyGPlates](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation).
- [PlateModelManager](https://pypi.org/project/plate-model-manager/) (note that installing [GPlately](https://gplates.github.io/gplately/stable/sphinx/html/index.html) also installs this).
- SciPy.
- And, on Windows platforms, optionally install [psutil](https://pypi.org/project/psutil/) so that this workflow can use CPU cores in the *background* (ie, below-normal priority).

## Usage

To generate the carbonate thickness grids you can either:

- load the Jupyter notebook `carbonate_sediment_thickness.ipynb` and run all cells, or
- type `python run_carbonate_sediment_thickness.py` in a console/terminal window.

In either case there are a bunch of top-level parameters that you can change/configure.
By default `use_all_cpu_cores` is set to `True` to run on all CPU cores
(otherwise it takes too long; up to 25 hours at 0.5 degree resolution using just a single core).
Note that you can increase the `grid_spacing` parameter to reduce the running time.

> **Note:** If you choose the Jupyter notebook *and* you edit a parameter *outside* the notebook
(such as inside the imported module *carbonate_sediment_thickness*) then you'll need to restart the notebook kernel
after each modification (or insert `reload(carbonate_sediment_thickness)` after `import carbonate_sediment_thickness`).

The location of the age and bathymetry will need to be changed to point to your local grids.

The age grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF.zip .

The bathymetry grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Wright_etal_2020_ESR/Grids/Paleobathymetry_RHCW18/ .

You can either use the supplied topological model (in local directory `input_data/topology_model/2019_v2/`), or use ``PlateModelManager`` to provide one, or provide your own. If you're using the supplied model then you don't need to do anything. If you're using ``PlateModelManager`` then you just need to specify the model name (and set ``using_local_model`` to ``False``). If you're providing your own model then you'll need to list your rotation and topology files in the `rotation_filenames` and `topology_filenames` variables (in `carbonate_sediment_thickness.ipynb` or `run_carbonate_sediment_thickness.py`) - and note that you'll likely need to use *absolute* paths in your filenames (unlike the supplied model).


## Latitude- and time-dependent correction (`latitude_factor.py`)

The carbonate-rate formulation in `carbonate_sediment_thickness.py` has no
temperature, productivity or carbonate-saturation-state control — only the
depth of each cell relative to a globally uniform CCD curve. Without an
additional correction this over-predicts pelagic carbonate thickness on
young, shallow, high-latitude crust (most visibly in the Arctic). The
module `latitude_factor.py` provides a multiplier in [0, 1] that the
inner reconstruction loop applies to the per-step carbonate rate.

The correction has two components:

1. **Spatial half-cosine taper in absolute palaeolatitude.** Factor = 1
   equatorward of `lat_full` (default 60°), decreasing smoothly to 0 at
   `lat_zero` (default 80°). The thresholds are calibrated against
   modern sediment-trap CaCO₃ flux compilations (Honjo et al. 2008), the
   modern deep-sea CaCO₃ atlas (Archer 1996), and the absence of pelagic
   carbonate in the deep Arctic (Stein 2008; Jutterström & Anderson
   2016). The factor is evaluated at the *paleo*-latitude of each
   reconstructed point at the current reconstruction time — so crust
   that originated at low latitude still gets credit for the carbonate
   it accumulated before drifting poleward.

2. **Temporal half-cosine ramp in palaeotime.** The strong latitudinal
   cut-off is an icehouse phenomenon and is not applied to pre-Antarctic-
   glaciation greenhouse intervals. The correction is switched off
   entirely at `t ≥ 34 Ma` (EOT-1 / Oi-1 Antarctic glaciation onset;
   Coxall et al. 2005; Miller et al. 2008), fully applied at `t ≤ 23 Ma`
   (end of the Oligocene CO₂ drawdown reconstructed by the CenCO2PIP
   consortium, Hönisch et al. 2023), and smoothly ramped in between.
   Direct evidence for pelagic carbonate at ≳55° S palaeolatitude through
   the pre-EOT interval comes from Maud Rise ODP 689/690 (Maastrichtian–
   Eocene oozes at ~65° S; Thomas 1990; Huber 1991), the Mentelle Basin
   IODP Expedition 369 sites U1513–U1516 (continuous pelagic carbonate
   across OAE 2 and the PETM at 54–60° S; Huber et al. 2019), and the
   Campbell Plateau / Otway Basin records that localise the carbonate
   shutdown to the EOT (Gallagher et al. 2020).

The two thresholds for each component are exposed as module-level
constants so they can be varied in sensitivity tests
(`HIGH_LATITUDE_FULL_PRESERVATION_DEG`, `HIGH_LATITUDE_ZERO_ACCUMULATION_DEG`,
`GLACIATION_ONSET_MA`, `FULL_CORRECTION_MA`). Suggested sensitivity sweeps:

- spatial: (50°, 70°), (60°, 80°), (65°, 85°)
- temporal: (34, 33), (34, 23), (34, 16)

`latitude_factor.py` contains the full reference list with DOIs. See also
`CHANGELOG.md` for the divergence from the upstream Dutkiewicz et al.
(2019) workflow and `email_to_developer.md` for the cover note prepared
for submission to the upstream repository.

## Reference

Dutkiewicz, A., Müller, R.D., Cannon, J., Vaughan, S. and Zahirovic, S., 2019, Sequestration and subduction of deep-sea carbonate in the global ocean since the Early Cretaceous. Geology, 47(1), pp.91-94. DOI:  https://doi.org/10.1130/G45424.1
