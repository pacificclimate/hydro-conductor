#!/usr/bin/env python

""" This script orchestrates a coupled Variable Infiltration Capacity (VIC) and Regional Glacier Model (RGM) run. """

import argparse
import csv
import collections
from collections import OrderedDict
from datetime import date, timedelta
import os
import subprocess
import sys

import numpy as np
import h5py

from snbparams import SnbParams
from vegparams import VegParams
from vic_globals import get_global_parms, update_global_parms, write_global_parms_file

vic_full_path = '/home/mfischer/code/vic/vicNl'  # should this be a command line parameter?
rgm_full_path = '/home/mfischer/code/rgm/rgm' # ditto?
temp_files_path = '/home/mfischer/vic_dev/out/testing/temp_out_files/' # ditto?
# set it as default = os.env(tmp)

BARE_SOIL_ID = '19'

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
    parser.add_argument('--bare-soil-root', action="store", dest="open_ground_root_zone_file", type=str, default=None, help = 'file name and path of one-line text file containing 6 custom root parameters for the bare soil vegetation type / HRU (same format as a vegetation tile line in the vegetation parameters file).  Default: 0.10  1.00  0.10  0.00  0.10  0.00')
    parser.add_argument('--glacier-root', action="store", dest="glacier_root_zone_file", type=str, default=None, help = 'file name and path of one-line text file containing 6 custom root parameters for the glacier vegetation type / HRU (same format as a vegetation tile line in the vegetation parameters file).  Default: 0.10  1.00  0.10  0.00  0.10  0.00')
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
    open_ground_root_zone_file = options.open_ground_root_zone_file
    glacier_root_zone_file = options.glacier_root_zone_file
    band_size = options.band_size

    if not open_ground_root_zone_file:
        open_ground_root_zone_parms = '0.10  1.00  0.10  0.00  0.10  0.00'
    else:
        open_ground_root_zone = np.loadtxt(options.open_ground_root_zone_file)
    if not glacier_root_zone_file:
        glacier_root_zone_parms = '0.10  1.00  0.10  0.00  0.10  0.00'
    else:
        glacier_root_zone_parms = np.loadtxt(options.glacier_root_zone_file)

    return vic_global_file, rgm_params_file, surf_dem_in_file, bed_dem_file, pixel_cell_map_file, init_glacier_mask_file, output_trace_files, glacier_root_zone_file, band_size

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
    """Reads in a Vegetation Parameter File (VPF) and parses out VIC grid cell IDs,
       and creates and populates a VegParams object to store vegetation parameters
       for referencing, modifying, and writing back out to a new VPF for subsequent VIC iterations.
    """
    vp = VegParams(veg_parm_file)
    return vp, vp.cell_ids

def get_snb_parms(snb_parm_file, num_snow_bands, band_size):
    """ Reads in a Snow Band Parameter File (SNB) and creates and populates a SnbParams object
        which also keeps track of the addition or loss of elevation bands vis-a-vis glacier evolution.
        The SnbParams object stores snow band (elevation) parameters for referencing, modifying, 
        and writing back out to a new SNB for subsequent VIC iterations.
    """
    sp = SnbParams(snb_parm_file, num_snow_bands, band_size)
    return sp

def update_band_area_fracs(cell_ids, cell_areas, snb_parms, veg_parms, num_snow_bands, band_size, pixel_to_cell_map,
                      surf_dem, num_rows_dem, num_cols_dem, glacier_mask):
    """ Calculates and updates the area fractions of elevation bands within VIC cells, and
        area fraction of glacier and open ground within VIC cells (broken down by elevation band).
        Calls VegParams::update_tiles() to update all vegetation parameter tiles (or add/delete them as needed)
    """
    all_pixel_elevs = {} # temporarily store pixel elevations in bins by band so the median can be calculated
    band_areas = {} # temporary count of pixels, a proxy for area, within each band
    glacier_areas = {} # temporary count of pixels landing within the glacier mask, a proxy for glacier area
    for cell in cell_ids:
        all_pixel_elevs[cell] = [[] for x in range(num_snow_bands)]
        band_areas[cell] = [0 for x in range(num_snow_bands)]
        glacier_areas[cell] = [0 for x in range(num_snow_bands)]
    for row in range(num_rows_dem):
        for col in range(num_cols_dem):
            cell = pixel_to_cell_map[row][col][0] # get the VIC cell this pixel belongs to
            if cell != 'NA':
                # Use the RGM DEM output to update the pixel median elevation in the pixel_to_cell_map
                pixel_elev = surf_dem[row][col]
                pixel_to_cell_map[row][col][1] = pixel_elev
                for band_idx, band in enumerate(snb_parms.cells[cell].band_map):
                    band_found = False
                    if int(pixel_elev) in range(band, band + band_size):
                        band_found = True
                        # Gather all pixel_elev values to update the median_elev of this band later
                        all_pixel_elevs[cell][band_idx].append(pixel_elev)
                        band_areas[cell][band_idx] += 1
                        if glacier_mask[row][col]:
                            glacier_areas[cell][band_idx] += 1
                        break
                if not band_found: # we have to introduce a new elevation band starting at int(pixel_elev - pixel_elev % band_size)
                    new_band_idx = snb_parms.create_band(cell, pixel_elev)
                    all_pixel_elevs[cell][new_band_idx].append(pixel_elev)
                    band_areas[cell][new_band_idx] += 1
                    if glacier_mask[row][col]:
                        glacier_areas[cell][new_band_idx] += 1
                        break
    # Update all band median elevations for all cells, delete unused bands
    for cell in cell_ids:
        for band_idx, band in enumerate(snb_parms.cells[cell].band_map):
            snb_parms.cells[cell].median_elev[band_idx] = np.median(all_pixel_elevs[cell][band_idx])
            # if no entries exist for this band in all_pixel_elevs, delete the band
            if not all_pixel_elevs[cell][band_idx]:
                delete_band(band)
    # Update all band area fractions and glacier area fractions for all bands in all cells
    for cell in cell_ids:
        print('cell: {}'.format(cell))
        for band_idx, band in enumerate(snb_parms.cells[cell].band_map):
            # Update area fraction for this cell[band]
            new_band_area_frac = float(band_areas[cell][band_idx]) / cell_areas[cell]
            snb_parms.cells[cell].area_fracs[band_idx] = new_band_area_frac
            if snb_parms.cells[cell].area_fracs[band_idx] < 0:
                print('update_band_area_fracs(): Error: Calculated a negative band area fraction for cell {}, band {}. Exiting.\n'.format(cell, band))
                sys.exit(0)
            # This will be used to update the Vegetation Parameter File
            new_glacier_area_frac = float(glacier_areas[cell][band_idx]) / cell_areas[cell]
            # If the glacier HRU area fraction has changed for this band then we need to update all area fractions
            if new_glacier_area_frac != veg_parms.cells[cell][band_idx].area_frac_glacier:
                if new_glacier_area_frac < 0:
                    print('update_band_area_fracs(): Error: Calculated a negative glacier area fraction for cell {}, band {}. Exiting.\n'.format(cell, band))
                    sys.exit(0)
                
                # Calculate new non-glacier area fraction for this band
                new_non_glacier_area_frac = new_band_area_frac - new_glacier_area_frac
                # Calculate new_residual area fraction
                new_residual_area_frac = new_non_glacier_area_frac - veg_parms.cells[cell][band_idx].area_frac_non_glacier
                # Use old proportions of vegetated areas for scaling their area fractions in this iteration
                veg_scaling_divisor = veg_parms.cells[cell][band_idx].area_frac_non_glacier - veg_parms.cells[cell][band_idx].area_frac_open_ground
                # Calculate new open ground area fraction
                new_open_ground_area_frac = np.max([0, (veg_parms.cells[cell][band_idx].area_frac_open_ground + new_residual_area_frac)])
                # Calculate the change in vegetated area
                delta_area_vegetated = np.min([0, (veg_parms.cells[cell][band_idx].area_frac_open_ground + new_residual_area_frac)])
                # Update glacier area fraction for this band
                veg_parms.cells[cell][band_idx].area_frac_glacier = new_glacier_area_frac
                # Update new open ground area fraction for this band
                veg_parms.cells[cell][band_idx].area_frac_open_ground = new_open_ground_area_frac
                # Update non-glacier area fraction for this band
                veg_parms.cells[cell][band_idx].area_frac_non_glacier = new_non_glacier_area_frac              
                # Update all vegetation type tiles' area fractions for this band (and add/remove glacier / open ground tiles if needed)
                veg_parms.update_tiles(cell, band, delta_area_vegetated, veg_scaling_divisor)


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

# Main program
def main():
    print('\n\nVIC + RGM ... together at last!')

    # Parse command line parameters
    vic_global_file, rgm_params_file, surf_dem_in_file, bed_dem_file, \
        pixel_cell_map_file, init_glacier_mask_file, output_trace_files, \
        glacier_root_zone_file, band_size = parse_input_parms()

    # Get all initial VIC global parameters from the global parameter file
    global_parms = get_global_parms(vic_global_file)

    # Get entire time range of coupled VIC-RGM run from the initial VIC global file
    start_date = date(global_parms['STARTYEAR'], global_parms['STARTMONTH'], global_parms['STARTDATE'])
    end_date = date(global_parms['ENDYEAR'], global_parms['ENDMONTH'], global_parms['ENDDATE'])
    # Set the initial year for the coupled VIC-RGM simulation
    year = start_date.year
    # Get the date that glacier accumulation is to start (at the end of VIC "spin-up")
    glacier_accum_start_date = date(global_parms['GLACIER_ACCUM_START_YEAR'],
                                    global_parms['GLACIER_ACCUM_START_MONTH'],
                                    global_parms['GLACIER_ACCUM_START_DAY'])

    # Initial VIC output state filename prefix is determined by STATENAME in the global file
    state_filename_prefix = global_parms['STATENAME']
    # Numeric code indicating a glacier vegetation tile (HRU)
    GLACIER_ID = global_parms['GLACIER_ID']

    # Get VIC vegetation parameters and grid cell IDs from initial Vegetation Parameter File
    veg_parm_file = global_parms['VEGPARAM']
    veg_parms, cell_ids = get_veg_parms(veg_parm_file)

    # Get VIC snow/elevation band parameters from initial Snow Band File
    num_snow_bands, snb_file = global_parms['SNOW_BAND'].split()
    num_snow_bands = int(num_snow_bands)
    snb_parms = get_snb_parms(snb_file, num_snow_bands, band_size)

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
    # NOTE: the following needs to be updated to handle initial case + general case. See whiteboard May 14
    update_band_area_fracs(cell_ids, cell_areas, snb_parms, veg_parms, num_snow_bands, band_size, pixel_to_cell_map,
                      surf_dem_initial, num_rows_dem, num_cols_dem, glacier_mask)
    
    # Calculate the initial non-glacier area fractions for all bands in all cells
    area_frac_non_glacier = veg_parms.init_non_glacier_area_fracs(snb_parms)
    
    # Update the vegetation parameters vis-a-vis the application of the initial glacier mask, and write to new temporary file temp_vpf
    #update_veg_parms(cell_ids, veg_parms, area_frac_bands, area_frac_glacier, residual_area_fracs)
    veg_parms.update(snb_parms, old_area_frac_glacier, None)
    temp_vpf = temp_files_path + 'vpf_temp_' + str(year) + '.txt'
    veg_parms.save(temp_vpf)
    
    # Update snow band parameters vis-a-vis the application of the initial glacier mask, and write to new temporary file temp_snb
    update_snb_parms(snb_parms, area_frac_bands)
    temp_snb = temp_files_path + 'snb_temp_' + str(year) + '.txt'
    write_snb_parms_file(temp_snb, snb_parms, area_frac_bands)
    

    # Run the coupled VIC-RGM model for the time range specified in the VIC global parameters file
    while year < end_date.year:
        print('\nRunning year: {}'.format(year))

        # 1. Write / Update temporary Global Parameters File, temp_gpf
        temp_gpf = temp_files_path + 'gpf_temp_' + str(year) + '.txt'
        #update_global_parms(global_parms, temp_vpf, temp_snb, num_snow_bands, \
        #    start_date, end_date, init_state_file, state_date)
        if year == start_date.year:
            # set global parameters for VIC "spin-up", from start_date.year to the end of the first glacier accumulation
            init_state_file = None
            temp_end_date = date(glacier_accum_start_date.year+1, glacier_accum_start_date.month, glacier_accum_start_date.day) - timedelta(days=1)
            update_global_parms(global_parms, temp_vpf, temp_snb, num_snow_bands, \
                start_date, temp_end_date, init_state_file, temp_end_date)
            # Fast-forward year to what it will be when VIC finishes its spin-up
            year = temp_end_date.year
        else:
            temp_start_date = date(year, glacier_accum_start_date.month, glacier_accum_start_date.day)
            temp_end_date = date(year+1, glacier_accum_start_date.month, glacier_accum_start_date.day) - timedelta(days=1)
            # set/create INIT_STATE parm as the last written VIC state_file (parm does not exist in the first read-in of global_parms)
            init_state_file = state_filename_prefix + "_" + str(year-1) + str(temp_end_date.month) + str(temp_end_date.day)
            update_global_parms(global_parms, temp_vpf, temp_snb, num_snow_bands, \
                temp_start_date, temp_end_date, init_state_file, temp_end_date)

        write_global_parms_file(global_parms, temp_gpf)
        print('invoking VIC with global parameter file {}'.format(temp_gpf))

        # 2. Run VIC for a year.  This will save VIC model state at the end of the year, along with a Glacier Mass Balance (GMB) polynomial for each cell
        subprocess.check_call([vic_full_path, "-g", temp_gpf], shell=False, stderr=subprocess.STDOUT)

        # 3. Open VIC NetCDF state file and get the most recent GMB polynomial for each grid cell being modeled
        state_file = state_filename_prefix + "_" + str(temp_end_date.year) + str(temp_end_date.month) + str(temp_end_date.day)
        print('opening VIC state file {}'.format(state_file))
        state = h5py.File(state_file, 'r+')
        gmb_polys = get_mass_balance_polynomials(state, state_file, cell_ids)
            
        # 4. Translate mass balances using grid cell GMB polynomials and current veg_parm_file into a 2D RGM mass balance grid (MBG)
        mass_balance_grid = mass_balances_to_rgm_grid(gmb_polys, pixel_to_cell_map, num_rows_dem, num_cols_dem, cell_ids)
        # write Mass Balance Grid to ASCII file to direct the RGM to use as input
        mbg_file = temp_files_path + 'mass_balance_grid_' + str(year) + '.gsa'
        write_grid_to_gsa_file(mass_balance_grid, mbg_file, num_cols_dem, num_rows_dem, dem_xmin, dem_xmax, dem_ymin, dem_ymax)

        # 5. Run RGM for one year, passing MBG, BDEM, SDEM
        #subprocess.check_call([rgm_full_path, "-p", rgm_params_file, "-b", bed_dem_file, "-d", sdem_file, "-m", mbg_file, "-o", temp_files_path, "-s", "0", "-e", "0" ], shell=False, stderr=subprocess.STDOUT)
        subprocess.check_call([rgm_full_path, "-p", rgm_params_file, "-b", bed_dem_file, "-d", surf_dem_in_file, "-m", mbg_file, "-o", temp_files_path, "-s", "0", "-e", "0" ], shell=False, stderr=subprocess.STDOUT)
        # remove temporary files if not saving for offline inspection
        if not output_trace_files:
            os.remove(mbg_file)
            os.remove(rgm_surf_dem_file)

        # 6. Read in new Surface DEM file from RGM output
        rgm_surf_dem_out = np.loadtxt(rgm_surf_dem_out_file, skiprows=5)
        temp_surf_dem_file = temp_files_path + 'rgm_surf_dem_out_' + str(year) + '.gsa'
        os.rename(rgm_surf_dem_out_file, temp_surf_dem_file)
        # this will be fed back into RGM on next time step
        surf_dem_in_file = temp_surf_dem_file

        # 7. Update glacier mask
        glacier_mask = update_glacier_mask(rgm_surf_dem_out, bed_dem, num_rows_dem, num_cols_dem)
        if output_trace_files:
            glacier_mask_file = temp_files_path + 'glacier_mask_' + str(year) + '.gsa'
            write_grid_to_gsa_file(glacier_mask, glacier_mask_file)
        
        # 8. Update areas of each elevation band in each VIC grid cell, and calculate area fractions
        update_band_area_fracs(cell_ids, cell_areas, snb_parms, veg_parms, num_snow_bands, band_size, pixel_to_cell_map, rgm_surf_dem_out, num_rows_dem, num_cols_dem, glacier_mask)

        # 9. Update vegetation parameters and write to new temporary file temp_vpf
        #update_veg_parms(cell_ids, veg_parms, area_frac_bands, area_frac_glacier, residual_area_fracs)
        veg_parms.update(snb_parms, old_area_frac_glacier, new_area_frac_glacier, area_frac_non_glacier)
        # Update old_area_frac_glacier to new_area_frac_glacier for next iteration
        old_area_frac_glacier = new_area_frac_glacier
        temp_vpf = temp_files_path + 'vpf_temp_' + str(year) + '.txt'
        veg_parms.save(temp_vpf)

        # 10. Update snow band parameters and write to new temporary file temp_snb
        update_snb_parms(snb_parms, area_frac_bands)
        temp_snb = temp_files_path + 'snb_temp_' + str(year) + '.txt'
        write_snb_parms_file(temp_snb, snb_parms, area_frac_bands)

        # 11 Update HRUs in VIC state file 
            # don't forget to close the state file
        
        # Increment the year for the next loop iteration
        year += 1

# Main program invocation.  
if __name__ == '__main__':
    main()
