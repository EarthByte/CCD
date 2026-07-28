
"""
    Copyright (C) 2018 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


##############################################################################
# Generate carbonate sediment thickness grids from age and bathymetry grids. #
##############################################################################



from call_system_command import call_system_command
import csv
import math
import multiprocessing
import os.path
import pygplates
from scipy.interpolate import interp1d
import sys


# Required pygplates version.
# Need version 0.50 to be able to use:
#   https://www.gplates.org/docs/pygplates/generated/pygplates.resolvedtopologicalboundary#pygplates.ResolvedTopologicalBoundary.reconstruct_point
PYGPLATES_VERSION_REQUIRED = pygplates.Version(0, 50)
#    
# Check the imported pygplates version.
if pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
    raise RuntimeError('Using pygplates version {0} but version {1} or greater is required'.format(
            pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))


#
# Input Parameters.
#


# Average ocean floor sediment density (kg/m^3), porosity (%) and decay (metres).
AVERAGE_OCEAN_FLOOR_SEDIMENT_DENSITY = 2647.0
AVERAGE_OCEAN_FLOOR_SEDIMENT_POROSITY = 0.66
AVERAGE_OCEAN_FLOOR_SEDIMENT_DECAY = 1333.0

# Dissolution distance (in metres) above the CCD.
# This distance above the CCD is where carbonate dissolution is 100%.
# Below this distance dissolution decreases linearly to 0% at the CCD.
CCD_DISSOLUTION_DISTANCE = 300.0


#
# Functions to read/generate input data.
#


# Generate uniformly spaced lat/lon grid of points in the specified lat/lon range.
# Points use GMT grid registration (start and end exactly on the dateline when using global lat/lon range).
def generate_input_points_grid(
        grid_spacing_degrees,
        latitude_range,   # (min, max) tuple
        longitude_range): # (min, max) tuple
    
    min_lat, max_lat = latitude_range
    min_lon, max_lon = longitude_range
    
    if (grid_spacing_degrees == 0 or
        max_lat <= min_lat or
        max_lon <= min_lon):
        return
    
    input_points = []
    
    # Data points start *on* 'min_lat' and 'min_lon'.
    # If '(max_lat - min_lat)' is an integer multiple of grid spacing then final latitude also lands on 'max_lat'.
    # If '(max_lon - min_lon)' is an integer multiple of grid spacing then final longitude also lands on 'max_lon'.
    num_latitudes = int(math.floor((max_lat - min_lat) / grid_spacing_degrees)) + 1
    num_longitudes = int(math.floor((max_lon - min_lon) / grid_spacing_degrees)) + 1
    for lat_index in range(num_latitudes):
        lat = min_lat + lat_index * grid_spacing_degrees
        
        for lon_index in range(num_longitudes):
            lon = min_lon + lon_index * grid_spacing_degrees
            
            input_points.append((lon, lat))
    
    return input_points


def gmt_grdtrack(
        input,
        grid_filename):
    """
    Samples a grid file at the specified locations.
    
    'input' is a list of (longitude, latitude) tuples where latitude and longitude are in degrees.
    
    Returns a list of float values.
    """
    
    # Create a multiline string (one line per lon/lat row).
    location_data = ''.join(
            ' '.join(str(item) for item in row) + '\n' for row in input)

    # The command-line strings to execute GMT 'grdtrack'.
    grdtrack_command_line = ["gmt", "grdtrack",
        "-G{0}".format(grid_filename),
        # Geographic input/output coordinates...
        "-fg",
        # Use linear interpolation, and avoid anti-aliasing...
        "-nl+a+bg+t0.5"]
    
    # Call the system command.
    stdout_data = call_system_command(grdtrack_command_line, stdin=location_data, return_stdout=True)

    # Extract the sampled values.
    output_values = []
    for line in stdout_data.splitlines():
        # Each line returned by GMT grdtrack contains "longitude latitude grid_value".
        # Note that if GMT returns "NaN" then we'll return float('nan').
        output_value = float(line.split()[2])
        output_values.append(output_value)
    
    return output_values


def read_curve(curve_filename):
    
    x_array = []
    y_array = []
    xmin, xmax = None, None
    with open(curve_filename, 'r') as curve_file:
        # Detect the delimiter (eg, tab, space, comma).
        csv_dialect = csv.Sniffer().sniff(curve_file.read(1024))
        curve_file.seek(0)
        # Read using the detected delimiter.
        curve_reader = csv.reader(curve_file, csv_dialect)
        for row in curve_reader:
            # Skip comments.
            if (row[0].startswith("#") or
                row[0].startswith(">")):
                continue

            if len(row) < 2:
                raise ValueError('Curve file "{0}" does not have at least 2 columns at line {1}.'.format(
                                 curve_filename,
                                 curve_reader.line_num - 1))

            x = float(row[0])
            y = float(row[1])
            
            x_array.append(x)
            y_array.append(y)
            
            # Track the y values of the two endpoints (ie, min/max x).
            if xmin is None or x < xmin:
                xmin = x
                y_at_xmin = y
            if xmax is None or x > xmax:
                xmax = x
                y_at_xmax = y
    
    # Convert curve to a piecewise linear interpolating 1D function that accepts x and returns y.
    # Return the y value at the boundary on bounds error instead of raising ValueError...
    return interp1d(x_array, y_array, bounds_error=False, fill_value=(y_at_xmin, y_at_xmax))


#
# Functions to write output data.
#


def write_xyz_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def write_grid_file_from_xyz(
        grid_filename,
        xyz_filename,
        grid_spacing,
        latitude_range,
        longitude_range):
    
    # The command-line strings to execute GMT 'nearneighbor'.
    # For example "nearneighbor thickness.xy -R-180/180/-90/90 -I1 -N4 -S1d -Gthickness.nc".
    gmt_command_line = [
        "gmt",
        "nearneighbor",
        xyz_filename,
        "-N4",
        "-S{0}d".format(1.5 * grid_spacing),
        "-I{0}".format(grid_spacing),
        "-R{0}/{1}/{2}/{3}".format(longitude_range[0], longitude_range[1], latitude_range[0], latitude_range[1]),
        # Use GMT gridline registration since our input point grid has data points on the grid lines.
        # Gridline registration is the default so we don't need to force pixel registration...
        # "-r",  # Force pixel registration since data points are at centre of cells.
        "-G{0}".format(grid_filename)]
    call_system_command(gmt_command_line)


def write_data(
        data,
        output_filename_prefix,
        grid_spacing,
        latitude_range,
        longitude_range):
    
    # Write XYZ file.
    xyz_filename = '{0}.xy'.format(output_filename_prefix)
    write_xyz_file(xyz_filename, data)
    
    # Write grid file.
    grid_filename = '{0}.nc'.format(output_filename_prefix)
    write_grid_file_from_xyz(grid_filename, xyz_filename, grid_spacing, latitude_range, longitude_range)


#
# Function to compact sediment thickness.
#

# Decompact total sediment thickness assuming a single lithology of the specified surface porosity and exponential decay.
#
# The decompacted (ie, surface) thickness d(z) of a compacted layer of height t(z) at depth z is:
#
#   d(z) = t(z) * (1 - porosity(z)) / (1 - porosity(0))
#
# ...due to the fact that the grain volume doesn't change with depth:
#
#   (1 - porosity(z1)) * t(z1) = (1 - porosity(z2)) * t(z2)
#
# ...and substituting z1=z and z2=0, and using d(z) = t(0) gives:
#
#   (1 - porosity(z)) * t(z) = (1 - porosity(0)) * d(z)
#   d(z) = t(z) * (1 - porosity(z)) / (1 - porosity(0))
#
# Where the porosity is:
#
#   porosity(z) = porosity(0) * exp(-c * z)
#
# The total decompacted thickness is:
#
#   D = Integral(d(z), z = 0 -> T)
#
# ...where T is the total compacted thickness.
#
#   D = Integral(dz * (1 - porosity(z)) / (1 - porosity(0)), z = 0 -> T)
#   D = Integral(dz * (1 - porosity(0) * exp(-c * z)) / (1 - porosity(0)), z = 0 -> T)
#
#     = (T + (porosity(0) / c) * (exp(-c * T) - 1)) / (1 - porosity(0))
#
# Rearranging gives...
#
#   T = D * (1 - porosity(0)) - (porosity(0) / c) * (exp(-c * T) - 1)
#     = a * exp(-c * T) + b
#
# ...where a and b are:
#
#    a = -porosity(0) / c
#
#    b = D * (1 - porosity(0)) - a
#
# So the compacted thickness T:
#
#    T = a * exp(-c * T) + b
#
# ...can be solved iteratively by repeatedly substituting the left hand side back into the right hand side
# until T converges on a solution. The initial T is chosen to be D (the decompacted thickness).
#
def compact_sediment_thickness(thickness, surface_porosity, porosity_exp_decay):
    
    # Constants 'a' and 'b' are calculated outside the iteration loop for efficiency.
    a = -surface_porosity / porosity_exp_decay
    b = thickness * (1 - surface_porosity) - a

    # Start out with initial estimate - choose the decompacted thickness.
    compacted_thickness = thickness

    # Limit the number of iterations in case we never converge.
    # Although should converge within around 20 iterations (for 1e-6 accuracy).
    for iteration in range(1000):
        new_compacted_thickness = a * math.exp(-porosity_exp_decay * compacted_thickness) + b
        
        # If we've converged within a tolerance then we're done.
        if math.fabs(new_compacted_thickness - compacted_thickness) < 1e-6:
            compacted_thickness = new_compacted_thickness
            break
        
        # The new thickness becomes the old thickness for the next loop iteration.
        compacted_thickness = new_compacted_thickness
    
    return compacted_thickness


def calc_carbonate_decompacted_sediment_thickness(
        input_points,
        age_grid_filenames_format,
        age_grid_anchor_plate_id,
        bathymetry_grid_filenames_format,
        bathymetry_grid_anchor_plate_id,
        bathymetry_filename_oldest_time,
        topological_model,
        ccd_curve,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve,
        carbonate_anchor_plate_id,
        time):
    """
    Calculates carbonate decompacted sediment thickness.
    
    Returns a list of tuples (lon, lat, age, bathymetry, carbonate_decompacted_sediment_thickness),
    Each tuple corresponds to an input point that has a valid age and a valid bathymetry.
    """

    # Convert points from (lon, lat) to pygplates.PointOnSphere.
    points = [pygplates.PointOnSphere(lat, lon) for lon, lat in input_points]

    # The rotation model is in the carbonate reference frame since that's the reference frame of the input points.
    rotation_model = topological_model.get_rotation_model()

    # Convert the points from the carbonate reference frame (default) to the age grid and bathymetry reference frames so
    # that we can sample the age and bathymetry grids at 'time' (using their respective reference frames).
    #
    # First define the reference frame conversion rotation (from carbonate reference frame C to bathymetry reference frame B).
    # Here 'P' is the plate ID of a point at time 't'.
    #
    # R = R(B->P, 0->t) * R(C->P, 0->t)^-1
    #   = R(B->701, 0->t) * R(701->P, 0->t) * [R(C->701, 0->t) * R(701->P, 0->t)]^-1
    #   = R(B->701, 0->t) * R(701->P, 0->t) * R(701->P, 0->t)^-1 * R(C->701, 0->t)^-1
    #   = R(B->701, 0->t) * R(C->701, 0->t)^-1
    #   = R(B->701, 0->t) * R(C->701, t->0)
    #
    # ...which is independent of 'P' and hence the same reference frame conversion rotation can be applied to all points.
    #
    # Arbitrary plate ID for converting between reference frames.
    non_reference_plate_id = 701
    # Get the rotation from 'time' to present day using the carbonate anchor plate ID (default anchor plate).
    rotation_from_carbonate_grid_reference_frame = rotation_model.get_rotation(0, non_reference_plate_id, time)
    # Get the rotation from present day back to 'time' using anchor plate 'age_grid_anchor_plate_id'.
    rotation_to_age_grid_reference_frame = rotation_model.get_rotation(
        time, non_reference_plate_id, 0, anchor_plate_id=age_grid_anchor_plate_id)
    # Get the rotation from present day back to 'time' using anchor plate 'bathymetry_grid_anchor_plate_id'.
    rotation_to_bathymetry_grid_reference_frame = rotation_model.get_rotation(
        time, non_reference_plate_id, 0, anchor_plate_id=bathymetry_grid_anchor_plate_id)
    #
    input_points_in_age_grid_reference_frame = []
    input_points_in_bathymetry_grid_reference_frame = []
    for point in points:
        input_lat_lon_point_in_age_grid_reference_frame = (
            rotation_to_age_grid_reference_frame * rotation_from_carbonate_grid_reference_frame * point)
        input_points_in_age_grid_reference_frame.append(
            input_lat_lon_point_in_age_grid_reference_frame.to_lat_lon()[::-1])  # swap (lat, lon)
        
        input_lat_lon_point_in_bathymetry_grid_reference_frame = (
            rotation_to_bathymetry_grid_reference_frame * rotation_from_carbonate_grid_reference_frame * point)
        input_points_in_bathymetry_grid_reference_frame.append(
            input_lat_lon_point_in_bathymetry_grid_reference_frame.to_lat_lon()[::-1])  # swap (lat, lon)
    
    # Age, distance and bathymetry filenames for the current time.
    age_grid_filename = age_grid_filenames_format.format(time)
    bathymetry_filename = bathymetry_grid_filenames_format.format(time)
    
    # Get the ages and bathymetries at the input point locations (in their respective reference frames).
    ages = gmt_grdtrack(input_points_in_age_grid_reference_frame, age_grid_filename)
    bathymetries = gmt_grdtrack(input_points_in_bathymetry_grid_reference_frame, bathymetry_filename)

    # Initialise all carbonate deposition to zero.
    # We will accumulate carbonate deposition only above the CCD.
    carbonate_decompacted_sediment_thicknesses = [0.0] * len(input_points)

    # Only reference points where there are both age *and* bathymetry values.
    point_indices = [point_index
                     for point_index in range(len(input_points))
                     if not math.isnan(ages[point_index]) and not math.isnan(bathymetries[point_index])]

    # Return early if no points have both age *and* bathymetry values.
    if len(point_indices) == 0:
        return

    # Time interval to reconstruct ocean points (from 'time' to each ocean point's time of appearance).
    # Currently set to 1 Myr.
    time_interval = 1.0
    
    reconstruction_time = time + time_interval
    reconstructed_ages = [(age - time_interval) for age in ages]
    reconstructed_point_indices = point_indices[:]
    reconstructed_points = [points[point_index] for point_index in point_indices]
    while reconstruction_time <= bathymetry_filename_oldest_time:
        # print(reconstruction_time); sys.stdout.flush()

        # Remove any ocean points that don't exist at the current reconstruction time
        # (ie, they appeared after current reconstruction time).
        index = 0
        while index < len(reconstructed_point_indices):
            reconstructed_point_index = reconstructed_point_indices[index]
            if reconstructed_ages[reconstructed_point_index] <= 0:
                del reconstructed_point_indices[index]
                del reconstructed_points[index]
                continue
            index += 1
        
        # Finished if all ocean points appeared after current reconstruction time.
        if not reconstructed_point_indices:
            break

        # Convert from the carbonate reference frame (default) to the bathymetry reference frame so that
        # we can sample the bathymetry grid at 'reconstruction_time' (using its reference frame).
        #
        # Get the rotation from 'reconstruction_time' to present day using the carbonate anchor plate ID (default anchor plate).
        rotation_from_carbonate_grid_reference_frame = rotation_model.get_rotation(0, non_reference_plate_id, reconstruction_time)
        # Get the rotation from present day back to 'reconstruction_time' using anchor plate 'bathymetry_grid_anchor_plate_id'.
        rotation_to_bathymetry_grid_reference_frame = rotation_model.get_rotation(
            reconstruction_time, non_reference_plate_id, 0, anchor_plate_id=bathymetry_grid_anchor_plate_id)

        # Find the topological plate/network containing each reconstructed point at 'reconstruction_time - time_interval'.
        reconstructed_point_locations = topological_model.topological_snapshot(
            reconstruction_time - time_interval).get_point_locations(reconstructed_points)
        
        # Reconstruct each point from 'reconstruction_time - time_interval' to 'reconstruction_time'.
        reconstructed_points_in_bathymetry_grid_reference_frame = []
        for index, reconstructed_point in enumerate(reconstructed_points):
            reconstructed_point_location = reconstructed_point_locations[index]

            # Reconstruct the current point to 'reconstruction_time' (from 'reconstruction_time - time_interval').
            if reconstructed_point_location.located_in_resolved_boundary():  # if point is inside a resolved boundary
                reconstructed_point = reconstructed_point_location.located_in_resolved_boundary().reconstruct_point(reconstructed_point, reconstruction_time)
            elif reconstructed_point_location.located_in_resolved_network():  # if point is inside a resolved network
                reconstructed_point = reconstructed_point_location.located_in_resolved_network().reconstruct_point(reconstructed_point, reconstruction_time)
            # else point is not in any resolved topologies, so leave it where it is (ie, don't reconstruct it over the current time interval).

            # Store the reconstructed point (so we can later continue reconstruction over the next time interval).
            reconstructed_points[index] = reconstructed_point
            
            # Convert to the bathymetry reference frame.
            reconstructed_point_in_bathymetry_grid_reference_frame = (
                rotation_to_bathymetry_grid_reference_frame * rotation_from_carbonate_grid_reference_frame * reconstructed_point)
            reconstructed_points_in_bathymetry_grid_reference_frame.append(
                reconstructed_point_in_bathymetry_grid_reference_frame.to_lat_lon()[::-1])  # swap (lat, lon)
        
        # Reconstructed bathymetry filename for the current reconstruction time.
        reconstructed_bathymetry_filename = bathymetry_grid_filenames_format.format(reconstruction_time)

        # Sample reconstructed bathymetry at the reconstructed locations.
        reconstructed_bathymetrys = gmt_grdtrack(reconstructed_points_in_bathymetry_grid_reference_frame, reconstructed_bathymetry_filename)
        
        # The CCD depth associated with the reconstructed time.
        ccd_depth_at_reconstruction_time = ccd_curve(reconstruction_time)

        # Reconstruct bathymetry for all ocean points that exist at the current reconstruction time.
        for index, reconstructed_point_index in enumerate(reconstructed_point_indices):
            
            # Sample bathymetry at current reconstruction time.
            #
            # Note the use of 'index' instead of 'reconstructed_point_index'
            bathymetry_at_reconstruction_time = reconstructed_bathymetrys[index]

            # print('  ', time, reconstruction_time, bathymetry_at_reconstruction_time); sys.stdout.flush()
            
            # If modelled bathymetry is above the CCD depth then add a 1 Myr time interval to the time spent above CDD.
            if (not math.isnan(bathymetry_at_reconstruction_time) and
                bathymetry_at_reconstruction_time > ccd_depth_at_reconstruction_time):
                
                # Linearly interpolate the carbonate sedimentation rate with depth between the maximum rate
                # at CCD_DISSOLUTION_DISTANCE metres above the CDD and zero (at the CCD).
                #
                # Because the dissolution of carbonate in the ocean is pressure and temperature dependent,
                # more and more carbonate gets dissolved as you go down. This effect ultimately causes the
                # CCD to occur, but is actually reduces the sedimentation rates as we go down from the MOR
                # towards deeper seafloor, and this phenomenon can be clearly seen in the sedimentation rates
                # in the graph computed in "Predicting sediment thickness on vanished ocean crust since 200 Ma".
                # You can see (in the previously mentioned sed rate graph) that the sed rate declines substantially
                # from the MOR depth at 0 Ma to the depth of the recent CCD (around 45-50 Ma) where the sed rate
                # flattens because below the CCD sed rate is no longer depth/age dependent. From published papers
                # it appears that it is a good assumption that the sed rate decreases linearly between the MOR depth
                # and the CCD. We know that it is zero for carbonates at the CCD, and the (maximum) value for the
                # MOR is obtained from a time-dependent carbonate sed rate file.
                # Of course if we have an oceanic plateau that is substantially shallower than the normal MOR this
                # rate might be even higher than what the graph implies for the MOR but we cannot do anything about that.
                # So what we want to do is make the sed rate linearly dependent on the depth range from the MOR to the CCD.
                #
                # UPDATE: It was found that starting at the MOR ridge (2.5km) for 100% dissolution was too conservative.
                #         Instead the dissolution zone is now always a fixed distance 'CCD_DISSOLUTION_DISTANCE' above the CCD.
                #
                # Note that previously we used the sed rate far away from cont margins, ie at 3km distance,
                # as a proxy for carbonate sed rate. This used the (maximum) value for the MOR from the
                # polynomial relationship for sed rate, ie we just plugged in 0 Ma and 3000 km. This assumed
                # the carbonate sedimentation rate equals the total sedimentation rate because the maximum mean
                # distance to margins (mentioned above) essentially removes sediment contributions from continents.
                # However now we have a fixed carbonate sed rate curve loaded from a file instead.
                #
                # Note that this is a decompacted thickness rate (ie, deposited at surface where no compaction takes place,
                # but porosity effects still apply, albeit only the *surface* porosity).
                #
                # NOTE: The factor of 10 converts cm/Ky to m/My (the file has units of cm/Ky).
                max_carbonate_decompacted_sediment_rate = 10 * max_carbonate_decomp_sed_rate_cm_per_ky_curve(reconstruction_time)
                carbonate_decompacted_sediment_rate = (
                    max_carbonate_decompacted_sediment_rate *
                    min(bathymetry_at_reconstruction_time - ccd_depth_at_reconstruction_time, CCD_DISSOLUTION_DISTANCE) /
                    CCD_DISSOLUTION_DISTANCE
                )
                carbonate_decompacted_sediment_thicknesses[reconstructed_point_index] += time_interval * carbonate_decompacted_sediment_rate

            # Decrease the age of crust as we step backward in time by one time interval.
            reconstructed_ages[reconstructed_point_index] -= time_interval

        # Increase the reconstruction time as we step backward in time by one time interval.
        reconstruction_time += time_interval

    # Combine carbonate decompacted sediment thickness with the lons, lats, ages and bathymetries.
    # And only include those points that have both age *and* bathymetry values (at 'time').
    return [input_points[point_index] + (ages[point_index], bathymetries[point_index], carbonate_decompacted_sediment_thicknesses[point_index])
            for point_index in point_indices]


def calc_sedimentation(
        input_points,  # List of (lon, lat) tuples,
        age_grid_filenames_format,
        age_grid_anchor_plate_id,
        bathymetry_grid_filenames_format,
        bathymetry_grid_anchor_plate_id,
        bathymetry_filename_oldest_time,
        topological_model,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        carbonate_anchor_plate_id,
        time):
    """
    Calculates decompacted and compacted carbonate thickness for each ocean basin grid point.
    Also determines when carbonate is currently being deposited (ie, only at the current time).
    
    Returns: A 3-tuple of lists each containing 3-tuples of (lon, lat, value).
             The values in the first list contain decompacted thickness.
             The values in the second list contain compacted thickness.
             The values in the third list contain a mask (value 1.0) for regions where
             carbonate is currently being deposited (ie, only at the current time).
    """
    
    # Read the CCD curve from the CCD file.
    # Returned curve is a function that accepts time and return depth.
    ccd_curve = read_curve(ccd_curve_filename)
    
    # Read the maximum carbonate decompacted sedimentation rate curve from the file.
    # Returned curve is a function that accepts time and return sedimentation rate (in cm/ky).
    max_carbonate_decomp_sed_rate_cm_per_ky_curve = read_curve(max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename)
    
    # Predict *decompacted* thickness.
    # We do this for all ocean basin points together since it's more efficient to sample reconstructed bathymetry for all points at once.
    lon_lat_age_bathymetry_carbonate_decompacted_sediment_thickness_list = calc_carbonate_decompacted_sediment_thickness(
        input_points,
        age_grid_filenames_format,
        age_grid_anchor_plate_id,
        bathymetry_grid_filenames_format,
        bathymetry_grid_anchor_plate_id,
        bathymetry_filename_oldest_time,
        topological_model,
        ccd_curve, max_carbonate_decomp_sed_rate_cm_per_ky_curve,
        carbonate_anchor_plate_id,
        time)
    if not lon_lat_age_bathymetry_carbonate_decompacted_sediment_thickness_list:
        return
    
    # The CCD at the current time.
    ccd_at_current_time = ccd_curve(time)
    
    # For each ocean basin point also generate *compacted* thickness.
    lon_lat_carbonate_decompacted_sediment_thickness_list = []
    lon_lat_carbonate_compacted_sediment_thickness_list = []
    lon_lat_carbonate_deposition_mask_list = []
    for lon, lat, age, bathymetry, carbonate_decompacted_sediment_thickness in lon_lat_age_bathymetry_carbonate_decompacted_sediment_thickness_list:
        
        # Compact thickness.
        carbonate_compacted_sediment_thickness = compact_sediment_thickness(
            carbonate_decompacted_sediment_thickness,
            AVERAGE_OCEAN_FLOOR_SEDIMENT_POROSITY,
            AVERAGE_OCEAN_FLOOR_SEDIMENT_DECAY)
        
        lon_lat_carbonate_decompacted_sediment_thickness_list.append((lon, lat, carbonate_decompacted_sediment_thickness))
        lon_lat_carbonate_compacted_sediment_thickness_list.append((lon, lat, carbonate_compacted_sediment_thickness))
        
        # See if carbonate sediment is being generated at the *current* time (if bathymetry above CCD).
        # This mask determines where on the current ocean floor carbonate is currently being deposited and
        # hence does not include the deposition history over the lifetime of the current parcel of ocean floor.
        if bathymetry > ccd_at_current_time:
            carbonate_deposition_mask = 1.0
            lon_lat_carbonate_deposition_mask_list.append((lon, lat, carbonate_deposition_mask))
    
    return (lon_lat_carbonate_decompacted_sediment_thickness_list,
            lon_lat_carbonate_compacted_sediment_thickness_list,
            lon_lat_carbonate_deposition_mask_list)


def calc_sedimentation_and_write_data(
        input_points,
        time,
        latitude_range,  # (min, max) tuple
        longitude_range, # (min, max) tuple
        grid_spacing,
        topological_model,
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
        carbonate_anchor_plate_id):
    
    print('Time: ', time)

    # Allow 'topological_model' to be passed as a tuple of (topology_filenames, rotation_filenames, anchor_plate_id).
    if not isinstance(topological_model, pygplates.TopologicalModel):
        topology_filenames, rotation_filenames, carbonate_anchor_plate_id = topological_model
        topological_model = pygplates.TopologicalModel(topology_filenames, rotation_filenames, anchor_plate_id=carbonate_anchor_plate_id)
    
    # Calculate carbonate decompacted and compacted sediment thickness at each input point that is in
    # the age and bathmetry grids (in unmasked regions of both grids).
    sediment_thickness_data = calc_sedimentation(
        input_points,
        age_grid_filenames_format,
        age_grid_anchor_plate_id,
        bathymetry_grid_filenames_format,
        bathymetry_grid_anchor_plate_id,
        bathymetry_filename_oldest_time,
        topological_model,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        carbonate_anchor_plate_id,
        time)
    if sediment_thickness_data is None:
        print('Ignoring time "{0}" - no grid points inside age and bathymetry grids.'.format(time), file=sys.stderr)
        return
    (carbonate_decompacted_sediment_thickness_data,
     carbonate_compacted_sediment_thickness_data,
     carbonate_deposition_mask_data) = sediment_thickness_data

    # Carbonate decompacted and compacted sediment thickness filenames without the filename extension.
    carbonate_decompacted_sediment_thickness_base_filename = (
        carbonate_decompacted_sediment_thickness_filename_prefix + '_{0}_{1}'.format(grid_spacing, time))
    carbonate_compacted_sediment_thickness_base_filename = (
        carbonate_compacted_sediment_thickness_filename_prefix + '_{0}_{1}'.format(grid_spacing, time))
    carbonate_deposition_mask_base_filename = (
        carbonate_deposition_mask_filename_prefix + '_{0}_{1}'.format(grid_spacing, time))
    
    # Output Carbonate decompacted and compacted sediment thicknesses to '.xy' files and '.nc' files.
    write_data(
        carbonate_decompacted_sediment_thickness_data,
        carbonate_decompacted_sediment_thickness_base_filename,
        grid_spacing, latitude_range, longitude_range)
    write_data(
        carbonate_compacted_sediment_thickness_data,
        carbonate_compacted_sediment_thickness_base_filename,
        grid_spacing, latitude_range, longitude_range)
    write_data(
        carbonate_deposition_mask_data,
        carbonate_deposition_mask_base_filename,
        grid_spacing, latitude_range, longitude_range)


#
# Parallel optimization
#


# Wraps around 'calc_sedimentation_and_write_data()' so can be used by multiprocessing.Pool.map()
# which requires a single-argument function.
def calc_sedimentation_and_write_data_parallel_pool_function(args):
    try:
        return calc_sedimentation_and_write_data(*args)
    except KeyboardInterrupt:
        pass


def low_priority():
    """ Set the priority of the process to below-normal."""

    import sys
    try:
        sys.getwindowsversion()
    except AttributeError:
        isWindows = False
    else:
        isWindows = True

    if isWindows:
        try:
            import psutil
        except ImportError:
            pass
        else:
            p = psutil.Process()
            p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    else:
        import os

        os.nice(1)


def calc_sedimentation_and_write_data_for_times(
        times,
        latitude_range,  # (min, max) tuple
        longitude_range, # (min, max) tuple
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
        use_all_cpu_cores):
    
    print('Starting...')
    
    # Generate a uniform global grid of lat/lon points.
    input_points = generate_input_points_grid(grid_spacing, latitude_range, longitude_range)
    
    # Either distribute each time iteration (over array of times) across all CPU cores, or
    # run time iteration loop serially (ie, use only one CPU core).
    if use_all_cpu_cores:

        # Instead of creating a topological model we pass in the data to create one.
        # This is because pyGPlates currently is a bit slow when pickling these objects.
        # Also, we can't reuse topological snapshots across different CPU cores since
        # a separate pygplates.TopologicalModel instance is created/pickled per CPU core.
        topological_model = topology_filenames, rotation_filenames, carbonate_anchor_plate_id
        
        # Split the workload across the CPUs.
        try:
            pool = multiprocessing.Pool(initializer=low_priority)
            pool_map_async_result = pool.map_async(
                calc_sedimentation_and_write_data_parallel_pool_function,
                (
                    (
                        input_points,
                        time,
                        latitude_range,
                        longitude_range,
                        grid_spacing,
                        topological_model,
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
                        carbonate_anchor_plate_id
                    ) for time in times
                ),
                1)  # chunksize

            # Apparently if we use pool.map_async instead of pool.map and then get the results
            # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
            # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
            try:
                pool_map_async_result.get(99999)
            except KeyboardInterrupt:
                # Note: 'finally' block below gets executed before returning.
                pass
        
        finally:
            pool.close()
            pool.join()

    else:  # process in serial...

        # Create a topological model up front to reuse topological snapshots across reconstruction times.
        topological_model = pygplates.TopologicalModel(
            topology_filenames, rotation_filenames, anchor_plate_id=carbonate_anchor_plate_id)

        # Iterate over times and generate grids for each time.
        for time in times:
            # Predict carbonate decompacted and compacted sediment thickness, and write results to output grids.
            calc_sedimentation_and_write_data(
                input_points,
                time,
                latitude_range,
                longitude_range,
                grid_spacing,
                topological_model,
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
                carbonate_anchor_plate_id)
    
    print('...finished.')
