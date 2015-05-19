#!/usr/bin/env python

""" This script orchestrates a coupled Variable Infiltration Capacity (VIC) and Regional Glacier Model (RGM) run. """

import argparse
import bisect
import csv
import collections
from collections import OrderedDict
from datetime import date, timedelta
import os
import subprocess
import sys
from warnings import warn

import numpy as np
import h5py
from dateutil.relativedelta import relativedelta

from conductor.vegparams import VegParams
from conductor.vic_globals import Global

vic_full_path = '/home/mfischer/code/vic/vicNl'  # should this be a command line parameter?
rgm_full_path = '/home/mfischer/code/rgm/rgm' # ditto?
temp_files_path = '/home/mfischer/vic_dev/out/testing/temp_out_files/' # ditto?
# set it as default = os.env(tmp)

BARE_SOIL_ID = '19'

one_year = relativedelta(years=+1)
one_day = relativedelta(days=+1)

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s]n' % message)
        self.print_help()
        sys.exit(2)

def parse_input_parms():
    # Get all global parameters 
    parser = MyParser()
    parser.add_argument('--g', action="store", dest="vic_global_file", type=str, help = 'file name and path of the VIC global parameters file')
    parser.add_argument('--rgm-params', action="store", dest="rgm_params_file", type=str, help = 'file name and path of the Regional Glacier Model (RGM) parameters file')
    parser.add_argument('--sdem', action="store", dest="surf_dem_file", type=str, help = 'file name and path of the initial Surface Digital Elevation Model (SDEM) file (GSA format)')
    parser.add_argument('--bdem', action="store", dest="bed_dem_file", type=str, help = 'file name and path of the Bed Digital Elevation Model (BDEM) file (GSA format)')
    parser.add_argument('--pixel-map', action="store", dest="pixel_cell_map_file", type=str, help = 'file name and path of the RGM Pixel to VIC Grid Cell mapping file')
    parser.add_argument('--glacier-mask', action="store", dest="init_glacier_mask_file", type=str, help = 'file name and path of the file containing the initial glacier mask (GSA format)')
    parser.add_argument('--trace-files', action="store_true", default=False, dest="trace_files", help = 'write out persistent GSA format surface DEMs and glacier masks, and 2D mass balance grid files, on each time step for offline inspection')
    parser.add_argument('--bare-soil-root', action="store", dest="bare_soil_root_parms_file", type=str, default=None, help = 'file name and path of one-line text file containing 6 custom root parameters for the bare soil vegetation type / HRU (same format as a vegetation tile line in the vegetation parameters file).  Default: 0.10  1.00  0.10  0.00  0.10  0.00')
    parser.add_argument('--glacier-root', action="store", dest="glacier_root_parms_file", type=str, default=None, help = 'file name and path of one-line text file containing 6 custom root parameters for the glacier vegetation type / HRU (same format as a vegetation tile line in the vegetation parameters file).  Default: 0.10  1.00  0.10  0.00  0.10  0.00')
    parser.add_argument('--band-size', action="store", dest="band_size", type=int, default=100, help = 'vertical size of VIC elevation bands in metres (default = 100m)')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    options = parser.parse_args()
    vic_global_file = options.vic_global_file
    rgm_params_file = options.rgm_params_file
    surf_dem_in_file = options.surf_dem_file
    bed_dem_file = options.bed_dem_file
    pixel_cell_map_file = options.pixel_cell_map_file
    init_glacier_mask_file = options.init_glacier_mask_file
    output_trace_files = options.trace_files
    bare_soil_root_parms_file = options.bare_soil_root_parms_file
    glacier_root_parms_file = options.glacier_root_parms_file
    band_size = options.band_size

    if not bare_soil_root_parms_file:
        bare_soil_root_parms = '0.10  1.00  0.10  0.00  0.10  0.00'
    else:
        bare_soil_root_parms = np.loadtxt(options.bare_soil_root_parms_file)
    if not glacier_root_parms_file:
        glacier_root_parms = '0.10  1.00  0.10  0.00  0.10  0.00'
    else:
        glacier_root_parms = np.loadtxt(glacier_root_parms_file)

    return vic_global_file, rgm_params_file, surf_dem_in_file, bed_dem_file, pixel_cell_map_file, init_glacier_mask_file, output_trace_files, glacier_root_parms_file, band_size

def read_gsa_headers(dem_file):
    """ Opens and reads the header metadata from a GSA Digital Elevation Map
        file, verifies agreement with the VIC-RGM mapping file metadata, and
        returns the x and y extents metadata
    """
    with open(dem_file, 'r') as f:
        # First line
        first_line = f.readline()
        assert first_line.startswith('DSAA'), 'read_gsa_headers({}): DSAA header on first line of DEM file was not found or is malformed.  DEM file does not conform to ASCII grid format.'.format(dem_file)
        # Second line
        num_cols, num_rows = f.readline().split()
        xmin, xmax = f.readline().split()
        ymin, ymax = f.readline().split()
        out_1 = [float(n) for n in (xmin, xmax, ymin, ymax)]
        out_2 = [int(x) for x in (num_rows, num_cols)]
    return out_1 + out_2

def get_rgm_pixel_mapping(pixel_map_file):
    """ Parses the RGM pixel to VIC grid cell mapping file and initialises a 2D
        grid of dimensions num_rows_dem x num_cols_dem (matching the RGM pixel
        grid), each element containing a list with the VIC cell ID associated
        with that RGM pixel and its median elevation
    """
    cell_areas = {}
    headers = {}
    
    with open(pixel_map_file, 'r') as f:

        # Read the number of columns and rows (order is unimportant)
        for _ in range(2):
            key, value = f.readline().split(None, 1)
            headers[key] = value

        num_cols_dem = int(headers['NCOLS'])
        num_rows_dem = int(headers['NROWS'])
        # create an empty two dimensional list
        pixel_to_cell_map = [[None] * num_cols_dem] * num_rows_dem 

        _ = f.readline() # Consume the column headers
        
        for line in f:
            # NOTE: column 3 is the elevation band; not used now, but maybe in future?
            _, row_num, col_num, _, median_elev, cell_id = line.split()
            pixel_to_cell_map[int(row_num)][int(col_num)] = (cell_id, median_elev)
            # Increment the pixel-granularity area within the grid cell
            if cell_id in cell_areas:
                cell_areas[cell_id] += 1
            else:
                cell_areas[cell_id] = 1

    return pixel_to_cell_map, num_rows_dem, num_cols_dem, cell_areas

def get_mass_balance_polynomials(state, state_file, cell_ids):
    """ Extracts the Glacier Mass Balance polynomial for each grid cell from an open VIC state file """
    gmb_info = state['GLAC_MASS_BALANCE_INFO'][0]
    cell_count = len(gmb_info)
    if cell_count != len(cell_ids):
        print('get_mass_balance_polynomials: The number of VIC cells ({}) read from the state file {} and those read from the vegetation parameter file ({}) disagree. Exiting.\n'.format(cell_count, state_file, len(cell_ids)))
        sys.exit(0)
    gmb_polys = {}
    for i in range(cell_count):
        cell_id = str(int(gmb_info[i][0]))
        if cell_id not in cell_ids:
            print('get_mass_balance_polynomials: Cell ID {} was not found in the list of VIC cell IDs read from the vegetation parameters file. Exiting.\n'.format(cell_id))
            sys.exit(0)
        gmb_polys[cell_id] = [gmb_info[i][1], gmb_info[i][2], gmb_info[i][3]]
    return gmb_polys

def mass_balances_to_rgm_grid(gmb_polys, pixel_to_cell_map, num_rows_dem, num_cols_dem, cell_ids):
    """ Translate mass balances from grid cell GMB polynomials to 2D RGM pixel grid to use as one of the inputs to RGM """
    mass_balance_grid = [[0 for x in range(num_cols_dem)] for x in range(num_rows_dem)]
    try:
        for row in range(num_rows_dem):
            for col in range(num_cols_dem):
                pixel = pixel_to_cell_map[row][col]
                # band = pixel[0]
                cell_id = pixel[0]
                median_elev = pixel[1] # read most recent median elevation of this pixel
                # only grab pixels that fall within a VIC cell
                if cell_id != 'NA':
                    # check that the cell_id agrees with what was read from the veg_parm_file
                    if cell_id not in cell_ids:
                        print('mass_balances_to_rgm_grid: Cell ID {} was not found in the list of VIC cell IDs read from the vegetation parameters file. Exiting.\n'.format(cell_id))
                        sys.exit(0)
                    mass_balance_grid[row][col] = gmb_polys[cell_id][0] + median_elev * (gmb_polys[cell_id][1] + median_elev * gmb_polys[cell_id][2])
    except:
        print('mass_balances_to_rgm_grid: Error while processing pixel {} (row {} column {})'.format(pixel, row, col))
    return mass_balance_grid

def write_grid_to_gsa_file(grid, outfilename, num_cols_dem, num_rows_dem, dem_xmin, dem_xmax, dem_ymin, dem_ymax):
    """ Writes a 2D grid to ASCII file in the input format expected by the RGM for DEM and mass balance grids """
    zmin = np.min(grid)
    zmax = np.max(grid)
    header_rows = [['DSAA'], [num_cols_dem, num_rows_dem], [dem_xmin, dem_xmax], [dem_ymin, dem_ymax], [zmin, zmax]]
    with open(outfilename, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ')
        for header_row in header_rows:
            writer.writerow(header_row)
        for row in grid:
            writer.writerow(row)

def get_veg_parms(veg_parm_file):
    """Reads in a Vegetation Parameter File and parses out VIC grid cell IDs,
       as well as an ordered dict of all vegetation parameters,
       grouped by elevation band index
    """
    vp = VegParams(veg_parm_file)
    return vp, vp.cell_ids

def update_band_areas(cell_ids, cell_areas, band_map, num_snow_bands, band_size, pixel_to_cell_map,
                      surf_dem, num_rows_dem, num_cols_dem, glacier_mask):
    """ Calculates the area fractions of elevation bands within VIC cells, and
        area fraction of glacier within VIC cells (broken down by elevation
        band)
    """
    band_areas = {}
    glacier_areas = {}
    area_frac_bands = {}
    area_frac_glacier = {}
    for cell in cell_ids:
        band_areas[cell] = [0 for x in range(num_snow_bands)]
        glacier_areas[cell] = [0 for x in range(num_snow_bands)]
        area_frac_bands[cell] = [0 for x in range(num_snow_bands)]
        area_frac_glacier[cell] = [0 for x in range(num_snow_bands)]
    for row in range(num_rows_dem):
        for col in range(num_cols_dem):
            cell = pixel_to_cell_map[row][col][0] # get the VIC cell this pixel belongs to
            if cell != 'NA':
                # Use the RGM DEM output to update the pixel median elevation in the pixel_to_cell_map
                pixel_elev = surf_dem[row][col]
                pixel_to_cell_map[row][col][1] = pixel_elev
                for band_idx, band in enumerate(band_map[cell]):
                    if int(pixel_elev) in range(band, band + band_size):
                        band_areas[cell][band_idx] += 1
                        if glacier_mask[row][col]:
                            glacier_areas[cell][band_idx] += 1
                        break

    print('Area fractions of elevation bands within VIC cells: {}'.format(band_areas))
    print('Area fraction of glacier within VIC cells, per elevation band: {}'.format(glacier_areas))
    print('\n')
    for cell in cell_ids:
        print('cell: {}'.format(cell))
        print('cell_areas[{}] = {}'.format(cell, cell_areas[cell]))
        print('\n')
        for band_idx, band in enumerate(band_map[cell]):
            print('band: {}, band_idx: {}'.format(band, band_idx))
            print('band_areas[{}][{}] = {}'.format(cell, band_idx, band_areas[cell][band_idx]))
            print('glacier_areas[{}][{}] = {}'.format(cell, band_idx, glacier_areas[cell][band_idx]))
            # This will be used to update the Snow Band File
            area_frac_bands[cell][band_idx] = float(band_areas[cell][band_idx]) / cell_areas[cell]
            if area_frac_bands[cell][band_idx] < 0:
                print('update_band_areas(): Error: Calculated a negative band area fraction for cell {}, band {}. Exiting.\n'.format(cell, band))
                sys.exit(0)
            # This will be used to update the Vegetation Parameter File
            area_frac_glacier[cell][band_idx] = float(glacier_areas[cell][band_idx]) / cell_areas[cell]
            if area_frac_glacier[cell][band_idx] < 0:
                print('update_band_areas(): Error: Calculated a negative glacier area fraction for cell {}, band {}. Exiting.\n'.format(cell, band))
                sys.exit(0)
            print('area_frac_bands[{}][{}] = {}'.format(cell, band_idx, area_frac_bands[cell][band_idx]))
            print('area_frac_glacier[{}][{}] = {}'.format(cell, band_idx, area_frac_glacier[cell][band_idx]))
            print('\n')
    return area_frac_bands, area_frac_glacier

def get_snb_parms(snb_file, num_snow_bands):
    """ Reads in a Snow Band File and outputs an ordered dict:
    {'cell_id_0' : [area_frac_band_0,...,area_frac_band_N],[median_elev_band_0,...,median_elev_band_N],[Pfactor_band_0,...,Pfactor_band_N]], 'cell_id_1' : ..."""
    snb_parms = OrderedDict()
    with open(snb_file, 'r') as f:
        for line in f:
            #print('snb file line: {}'.format(line))
            split_line = line.split()
            num_columns = len(split_line)
            cell_id = split_line[0]
            if num_columns != 3*num_snow_bands + 1:
                print('get_snb_parms(): Error: Number of columns ({}) in snow band file {} is incorrect for the given number of snow bands ({}) given in the global parameter file (should be 3 * num_snow_bands + 1). Exiting.\n'.format(num_columns, snb_file, num_snow_bands))
                sys.exit(0)
            snb_parms[cell_id] = [[float(x) for x in split_line[1 : num_snow_bands+1]],[int(x) for x in split_line[num_snow_bands+1 : 2*num_snow_bands+1]],[float(x) for x in split_line[2*num_snow_bands+1 : 3*num_snow_bands+1]]]
    return snb_parms

def create_band_map(snb_parms, band_size, cell_ids):
    """ Takes a dict of Snow Band parameters and identifies and creates a list
        of elevation bands for each grid cell of band_size width in meters
    """
    #TODO: create command line parameter to allow for extra bands to be specified as padding above/below existing ones, to allow for glacier growth(/slide?)
    band_map = {}
    for cell in cell_ids:
        band_map[cell] = [ int(band - band % band_size) for band in snb_parms[cell][1] ]
    return band_map

def update_glacier_mask(sdem, bdem, num_rows_dem, num_cols_dem):
    """ Takes output Surface DEM from RGM and uses element-wise differencing 
        with the Bed DEM to form an updated glacier mask 
    """
    diffs = sdem - bed_dem
    if np.any(diffs < 0):
        print('update_glacier_mask: Error: subtraction of Bed DEM from output Surface DEM of RGM produced one or more negative values.  Exiting.\n')
        sys.exit(0)
    glacier_mask = np.zeros((num_rows_dem, num_cols_dem))
    glacier_mask[diffs > 0] = 1
    return glacier_mask

def update_snb_parms(snb_parms, area_frac_bands):
    print('update_snb_parms...')
    for cell in snb_parms:
        #print('cell: {}'.format(cell)
        snb_parms[cell][0] = area_frac_bands[cell]
        print(' '.join(map(str, snb_parms[cell][0])))

def write_snb_parms_file(temp_snb, snb_parms, area_frac_bands):
    """ Writes current (updated) snow band parameters to a new temporary
        Snow Band File for feeding back into VIC
    """
    print('write_snb_parms_file...')
    #NOTE: the following does not work for the list of lists that is snb_parms
    # with open(temp_snb, 'w') as f:
    #     writer = csv.writer(f, delimiter=' ')
    #     for cell in snb_parms: 
    #         print('area_frac_bands[{}]: {}'.format(cell, area_frac_bands[cell]))
    #         print('snb_parms[{}][1:3]: {}'.format(cell, snb_parms[cell][1:3]))                  
    #         line = [ cell, area_frac_bands[cell] ] + snb_parms[cell][1:3]
    #         writer.writerow(line)
    with open(temp_snb, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        for cell in snb_parms:
            line = []
            line.append(cell)
            for area_frac in area_frac_bands[cell]:
                line.append(area_frac)
            for band_frac in snb_parms[cell][1]:
                line.append(band_frac) # append existing median elevations
            for pfactor in snb_parms[cell][2]:
                line.append(pfactor) # append existing Pfactor values
            writer.writerow(line)

def run_ranges(startdate, enddate, glacier_start):
    '''Generator which yields date ranges (a 2-tuple) that represent
       times at which to begin and end a VIC run.
       startdate and enddate are the overall boundaries of the
       simulation. glacier_start is a date some time after startdate
       that *should* be aligned with the water year (this code mostly
       assumes that that is the case).
       arguments (1950/01/01, 1995/12/31, 1955/10/01) would yield a
       sequence as follows:
      (1950/01/01, 1956/09/30)
      (1956/10/01, 1957/09/30)
      (1957/10/01, 1958/09/30)
      (1958/10/01, 1959/09/30)
      ...
      (1994/10/01, 1995/09/30)
      Note that the sequence will and should end before the specified
      end date, aligned with the water year.
    '''
    if not ((glacier_start.month, glacier_start.day) == (10, 1)):
        warn("run_ranges assumes that glacier_start is aligned to the water year"
             "(October 1 - September 30). You have configured glacier_start to"
             "{}. Only do this if you *really* know what you're doing"\
             .format(glacier_start.isoformat()))

    # First iteration doesn't include glaciers and ends at the water year
    t0 = startdate
    tn = glacier_start - one_day + one_year
    yield t0, tn
    # Secord itration aligned to water year, include glaciers
    # All subsequent iterations are aligned to water year and include glacier
    while tn < enddate:
        t0 = tn + one_day
        tn = t0 - one_day + one_year
        # don't let the simulation proceed beyond the enddate
        tn = min(tn, enddate)
        yield t0, tn

# Main program
def main():
    print('\n\nVIC + RGM ... together at last!')

    # Parse command line parameters
    vic_global_file, rgm_params_file, surf_dem_in_file, bed_dem_file, \
        pixel_cell_map_file, init_glacier_mask_file, output_trace_files, \
        glacier_root_parms_file, band_size = parse_input_parms()

    # Get all initial VIC global parameters from the global parameter file
    global_parms = Global(vic_global_file)

    assert global_parms.state_format == 'NETCDF', \
        "{} only supports NetCDF input statefile input as opposed "\
        "to the specified {}. Please change change this in your "\
        "global file {}".format(__name__, global_parms.state_format,
                                vic_global_file)

    # Initial VIC output state filename prefix is determined by STATENAME in the global file
    state_filename_prefix = global_parms.statename

    # Get VIC vegetation parameters and grid cell IDs from initial Vegetation Parameter File
    veg_parms, cell_ids = get_veg_parms(global_parms.vegparm)

    # Get VIC snow/elevation band parameters from initial Snow Band File
    num_snow_bands, snb_file = global_parms.snow_band.split()
    num_snow_bands = int(num_snow_bands)
    snb_parms = get_snb_parms(snb_file, num_snow_bands)

    # Get list of elevation bands for each VIC grid cell
    band_map = create_band_map(snb_parms, band_size, cell_ids)

    # The RGM will always output a DEM file of the same name (if running RGM for a single year at a time)
    rgm_surf_dem_out_file = temp_files_path + 's_out_00001.grd'

    # Open and read VIC-grid-to-RGM-pixel mapping file
    # pixel_to_cell_map is a list of dimensions num_rows_dem x num_cols_dem, each element containing a VIC grid cell ID
    pixel_to_cell_map, num_rows_dem, num_cols_dem, cell_areas = get_rgm_pixel_mapping(pixel_cell_map_file)

    # Get DEM xmin, xmax, ymin, ymax metadata of Bed DEM and check file header validity     
    dem_xmin, dem_xmax, dem_ymin, dem_ymax, num_rows, num_cols = read_gsa_headers(bed_dem_file)
    # Verify number of columns & rows agree with what's stated in the pixel_to_cell_map_file
    assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch of stated dimension(s) \
        between Bed DEM in {} (num rows: {}, num columns: {}) and the VIC-grid-to-RGM-pixel map in {} \
        (num rows: {}, num columns: {}). Exiting.\n'.format(bed_dem_file, num_rows, num_cols, pixel_cell_map_file, num_rows_dem, num_cols_dem)
    # Read in the provided Bed Digital Elevation Map (BDEM) file to 2D bed_dem array
    bed_dem = np.loadtxt(bed_dem_file, skiprows=5)

    # Check header validity of Surface DEM file
    _, _, _, _, num_rows, num_cols = read_gsa_headers(surf_dem_in_file)
    # Verify number of columns & rows agree with what's stated in the pixel_to_cell_map_file
    assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch of stated dimension(s) \
        between Surface DEM in {} (num rows: {}, num columns: {}) and the VIC-grid-to-RGM-pixel map in {} \
        (num rows: {}, num columns: {}). Exiting.\n'.format(surf_dem_in_file, num_rows, num_cols, pixel_cell_map_file, num_rows_dem, num_cols_dem)
    # Read in the provided Surface Digital Elevation Map (SDEM) file to 2D surf_dem array
    surf_dem_initial = np.loadtxt(surf_dem_in_file, skiprows=5)

    # Check header validity of initial Glacier Mask file
    _, _, _, _, num_rows, num_cols = read_gsa_headers(init_glacier_mask_file)
    # Verify number of columns & rows agree with what's stated in the pixel_to_cell_map_file
    assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch of stated dimension(s) \
        between Glacier Mask in {} (num rows: {}, num columns: {}) and the VIC-grid-to-RGM-pixel map in {} \
        (num rows: {}, num columns: {}). Exiting.\n'.format(init_glacier_mask_file, num_rows, num_cols, pixel_cell_map_file, num_rows_dem, num_cols_dem)
    # Read in the provided initial glacier mask file to 2D glacier_mask array 
    glacier_mask = np.loadtxt(init_glacier_mask_file, skiprows=5)

    # Apply the initial glacier mask and modify the band and glacier area fractions accordingly
    area_frac_bands, area_frac_glacier = update_band_areas(cell_ids, cell_areas, band_map, num_snow_bands, band_size, pixel_to_cell_map,
                      surf_dem_initial, num_rows_dem, num_cols_dem, glacier_mask)
    input('')
    # Calculate the initial residual (i.e. non-glacier) area fractions for all bands in all cells
    #residual_area_fracs = init_residual_area_fracs(cell_ids, veg_parms, snb_parms)
    veg_parms.init_residual_area_fracs(snb_parms)
    input('')
    # Update the vegetation parameters vis-a-vis the application of the initial glacier mask, and write to new temporary file temp_vpf
    #update_veg_parms(cell_ids, veg_parms, area_frac_bands, area_frac_glacier, residual_area_fracs)
    veg_parms.update(area_frac_bands, area_frac_glacier)
    temp_vpf = temp_files_path + 'vpf_temp_' + global_parms.startdate.isoformat() + '.txt'
    veg_parms.save(temp_vpf)
    input('')
    # Update snow band parameters vis-a-vis the application of the initial glacier mask, and write to new temporary file temp_snb
    update_snb_parms(snb_parms, area_frac_bands)
    temp_snb = temp_files_path + 'snb_temp_' + global_parms.startdate.isoformat() + '.txt'
    write_snb_parms_file(temp_snb, snb_parms, area_frac_bands)
    input('')

    # Run the coupled VIC-RGM model for the time range specified in the VIC global parameters file
    time_iterator = run_ranges(global_parms.startdate,
                               global_parms.enddate,
                               global_parms.glacier_start)
    for start, end in time_iterator:
        print('\nRunning VIC from {} to {}'.format(start, end))

        # 1. Write / Update temporary Global Parameters File, temp_gpf
        temp_gpf = temp_files_path + 'gpf_temp_{}.txt'.format(start.isoformat())

        # set global parameters for this VIC run
        global_parms.vegparm = temp_vpf
        global_parms.snow_band = '{} {}'.format(temp_snb, num_snow_bands)
        global_parms.startdate = start
        global_parms.enddate = end
        global_parms.statedate = end

        global_parms.write(temp_gpf)
        print('invoking VIC with global parameter file {}'.format(temp_gpf))

        # 2. Run VIC for a year.  This will save VIC model state at the end of the year, along with a Glacier Mass Balance (GMB) polynomial for each cell
        subprocess.check_call([vic_full_path, "-g", temp_gpf], shell=False, stderr=subprocess.STDOUT)

        # 3. Open VIC NetCDF state file and get the most recent GMB polynomial for each grid cell being modeled
        state_file = state_filename_prefix + "_" + start.isoformat() ## FIXME
        print('opening VIC state file {}'.format(state_file))
        state = h5py.File(state_file, 'r+')
        gmb_polys = get_mass_balance_polynomials(state, state_file, cell_ids)
            
        # 4. Translate mass balances using grid cell GMB polynomials and current veg_parm_file into a 2D RGM mass balance grid (MBG)
        mass_balance_grid = mass_balances_to_rgm_grid(gmb_polys, pixel_to_cell_map, num_rows_dem, num_cols_dem, cell_ids)
        # write Mass Balance Grid to ASCII file to direct the RGM to use as input
        mbg_file = temp_files_path + 'mass_balance_grid_' + start.isoformat() + '.gsa'
        write_grid_to_gsa_file(mass_balance_grid, mbg_file, num_cols_dem, num_rows_dem, dem_xmin, dem_xmax, dem_ymin, dem_ymax)

        # 5. Run RGM for one year, passing MBG, BDEM, SDEM
        subprocess.check_call([rgm_full_path, "-p", rgm_params_file, "-b", bed_dem_file, "-d", surf_dem_in_file, "-m", mbg_file, "-o", temp_files_path, "-s", "0", "-e", "0" ], shell=False, stderr=subprocess.STDOUT)
        # remove temporary files if not saving for offline inspection
        if not output_trace_files:
            os.remove(mbg_file)
            os.remove(rgm_surf_dem_file)

        # 6. Read in new Surface DEM file from RGM output
        rgm_surf_dem_out = np.loadtxt(rgm_surf_dem_out_file, skiprows=5)
        temp_surf_dem_file = temp_files_path + 'rgm_surf_dem_out_' + start.isoformat() + '.gsa'
        os.rename(rgm_surf_dem_out_file, temp_surf_dem_file)
        # this will be fed back into RGM on next time step
        surf_dem_in_file = temp_surf_dem_file

        # 7. Update glacier mask
        glacier_mask = update_glacier_mask(rgm_surf_dem_out, bed_dem, num_rows_dem, num_cols_dem)
        if output_trace_files:
            glacier_mask_file = temp_files_path + 'glacier_mask_' + start.isoformat() + '.gsa'
            write_grid_to_gsa_file(glacier_mask, glacier_mask_file)
        
        # 8. Update areas of each elevation band in each VIC grid cell, and calculate area fractions
        area_frac_bands, area_frac_glacier = update_band_areas(cell_ids, cell_areas, band_map, num_snow_bands, band_size, pixel_to_cell_map, rgm_surf_dem_out, num_rows_dem, num_cols_dem, glacier_mask)

        # 9. Update vegetation parameters and write to new temporary file temp_vpf
        #update_veg_parms(cell_ids, veg_parms, area_frac_bands, area_frac_glacier, residual_area_fracs)
        veg_parms.update(area_frac_bands, area_frac_glacier, residual_area_fracs)
        temp_vpf = temp_files_path + 'vpf_temp_' + start.isoformat() + '.txt'
        veg_parms.save(temp_vpf)

        # 10. Update snow band parameters and write to new temporary file temp_snb
        update_snb_parms(snb_parms, area_frac_bands)
        temp_snb = temp_files_path + 'snb_temp_' + start.isoformat() + '.txt'
        write_snb_parms_file(temp_snb, snb_parms, area_frac_bands)

        # 11 Update HRUs in VIC state file
            # don't forget to close the state file

        # Get ready for the next loop
        global_parms.init_state = "{}_{}.txt".format(global_parms.statename, end.strftime("%Y%m%d"))
        global_parms.statedate = end

# Main program invocation.
if __name__ == '__main__':
    main()
