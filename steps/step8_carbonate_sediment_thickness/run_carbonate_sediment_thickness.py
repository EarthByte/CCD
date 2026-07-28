#
# Generate carbonate sediment thickness grids from paleobathymetry grids.
#
# See "carbonate_sediment_thickness.py" for details of the algorithm.
#

import carbonate_sediment_thickness
import glob
import os.path
import os


#
# Input Parameters.
#

# Grid spacing of ocean basin points (in degrees).
grid_spacing = 0.5

# Regional latitude and longitude range.
# Keep latitudes in the range [-90, 90] and longitudes the in range [-180, 180].
# Set to [-90, 90] and [-180, 180] to cover the entire globe.
min_lat, max_lat = -90, 90
min_lon, max_lon = -180, 180

# Times to generate sediment thickness grids.
# Must also have age and bathymetry grids at these times.
times = list(range(0, 181))

# Whether to use all CPU cores (parallel) or just one (serial).
# Note: Each process is set to a low priority so as not to interfere with your regular tasks.
use_all_cpu_cores = True

# The topological model used to assign plate IDs to ocean crust at paleo times (including crust subducted at present day).
#
# Currently the only builtin model is provided in the sub-directory 'input_data/topology_model/2019_v2/'.
#
# However you can provide your own topological model either by listing the rotation and topology files or
# using plate-model-manager (which will download a plate model to the 'plate-model-repo' directory).
using_local_model = True
if using_local_model:
    rotation_filenames = glob.glob(os.path.join('input_data', 'topology_model', '2019_v2', '*.rot'))
    topology_filenames = glob.glob(os.path.join('input_data', 'topology_model', '2019_v2', '*.gpmlz'))
else:
    plate_model_name = "Muller2019"
    # Obtain rotation and topology files from Muller et al. (2019) via the plate model manager.
    #
    from plate_model_manager import PlateModelManager
    plate_model_manager = PlateModelManager()
    plate_model = plate_model_manager.get_model(plate_model_name, data_dir="plate-model-repo")
    rotation_filenames = plate_model.get_rotation_model()
    topology_filenames = plate_model.get_topologies()

# The reference frame (anchor plate ID) of the *output* carbonate grids.
carbonate_anchor_plate_id = 701

# CCD (calcite compensation depth) curve filename.
# This file maps time to CCD depth (negative).
ccd_curve_filename = 'input_data/Boss_Wilkinson_1991_global_CCD.txt'

# Maximum carbonate decompacted sediment rate curve filename (in cm/ky).
# This file maps time to the maximum carbonate rate (at mid-ocean ridge depth; reduces to zero at CCD).
max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename = 'input_data/sed_rate_v6.txt'

# Age grid files.
#
# The format string to generate age grid filenames (using the paleo times in 'times').
# Use a string section like "{:.1f}" to for the paleo time. The ".1f" part means use the paleo time to one decimal place
# (see Python\'s str.format() function) such that a time of 100 would be substituted as "100.0".
# This string section will get replaced with each time in turn (to generate the actual age grid filenames).
age_grid_filenames_format = 'C:/Users/jcann/Development/Usyd/data/EarthbytePlateModel/Muller_etal_2022/mantle-ref-frame-oceanic-crustal-agegrids_v1.2/Muller2022_SEAFLOOR_AGE_grid_{:.1f}Ma.nc'

# The reference frame (anchor plate ID) of the *input* age grids.
age_grid_anchor_plate_id = 0

# Bathymetry grid files.
#
# The format string to generate bathymetry grid filenames (using the paleo times in 'times').
# Use a string section like "{:.1f}" to for the paleo time. The ".1f" part means use the paleo time to one decimal place
# (see Python\'s str.format() function) such that a time of 100 would be substituted as "100.0".
# This string section will get replaced with each time in turn (to generate the actual bathymetry filenames).
bathymetry_grid_filenames_format = 'C:/Users/jcann/Development/Usyd/data/EarthbytePlateModel/Wright_etal_2020_ESR/Grids-Muller_etal_2019_v2.0/Paleobathymetry_RHCW18-Muller++_2019_v2.0/paleobathymetry_{:.0f}Ma.nc'

# The reference frame (anchor plate ID) of the *input* bathymetry grids.
bathymetry_grid_anchor_plate_id = 0

# Time of the oldest bathymetry grid.
bathymetry_filename_oldest_time = 230

# Location of output carbonate decompacted and compacted sediment thickness grids.
# The full path to thickness grids including the base filename
# (grid spacing, time and filename extension will get added later).
output_data_dir = 'sediment_thickness'
carbonate_decompacted_sediment_thickness_filename_prefix = os.path.join(output_data_dir, 'decompacted_sediment_thickness')
carbonate_compacted_sediment_thickness_filename_prefix = os.path.join(output_data_dir, 'compacted_sediment_thickness')
carbonate_deposition_mask_filename_prefix = os.path.join(output_data_dir, 'deposition_mask')


if __name__ == '__main__':
    
    # Create output directory if it doesn't exist.
    if not os.path.exists(output_data_dir):
        os.makedirs(output_data_dir)
    
    # All functionality is delegated to "carbonate_sediment_thickness.py" to enable parallel processing.
    carbonate_sediment_thickness.calc_sedimentation_and_write_data_for_times(
        times,
        (min_lat, max_lat),
        (min_lon, max_lon),
        grid_spacing,
        rotation_filenames,
        topology_filenames,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        age_grid_filenames_format,
        age_grid_anchor_plate_id,
        bathymetry_grid_filenames_format,
        bathymetry_grid_anchor_plate_id,
        bathymetry_filename_oldest_time,
        carbonate_decompacted_sediment_thickness_filename_prefix,
        carbonate_compacted_sediment_thickness_filename_prefix,
        carbonate_deposition_mask_filename_prefix,
        carbonate_anchor_plate_id,
        use_all_cpu_cores)
