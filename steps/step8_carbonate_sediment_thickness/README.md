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


## Reference

Dutkiewicz, A., Müller, R.D., Cannon, J., Vaughan, S. and Zahirovic, S., 2019, Sequestration and subduction of deep-sea carbonate in the global ocean since the Early Cretaceous. Geology, 47(1), pp.91-94. DOI:  https://doi.org/10.1130/G45424.1
