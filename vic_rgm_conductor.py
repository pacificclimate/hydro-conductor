#!/usr/bin/env python

""" This script orchestrates a coupled Variable Infiltration Capacity (VIC) and Regional Glacier Model (RGM) run. """

import argparse
import bisect
import csv
import collections
from collections import OrderedDict
from datetime import date
import os
import subprocess
import sys

import numpy as np
import h5py

import vegparams

vic_full_path = '/home/mfischer/code/vic/vicNl'  # should this be a command line parameter?
rgm_full_path = '/home/mfischer/code/rgm/rgm' # ditto?
temp_files_path = '/home/mfischer/vic_dev/out/testing/temp_out_files/' # ditto?
# set it as default = os.env(tmp)

BARE_SOIL_ID = '19'


# To have nested ordered defaultdicts
class OrderedDefaultdict(OrderedDict):
    # from: http://stackoverflow.com/questions/4126348/how-do-i-rewrite-this-function-to-implement-ordereddict/4127426#4127426
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or isinstance(args[0], collections.Callable)):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultdict, self).__init__(*args, **kwargs)
    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default
    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else ()
        return self.__class__, args, None, None, iter(self.items())

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
        pixel_to_cell_map = [[None] * num_col_dems] * num_rows_dem 

        _ = f.readline() # Consume the column headers
        
        for line in f:
            _, row_num, col_num, median_elev, cell_id = line.split()
            pixel_to_cell_map[row_num][col_num] = (cell_id, median_elev)
            # Increment the pixel-normalized area within the grid cell
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

def get_global_parms(global_parm_file):
    """ Parses the initial VIC global parameters file created by the user with the settings for the entire VIC-RGM run """
    global_parms = OrderedDefaultdict()
    n_outfile_lines = 0
    with open(global_parm_file, 'r') as f:
        for line in f:
            #print('line: {}'.format(line))
            if not line.isspace() and line[0] is not '#':
                split_line = line.split()
                #print('columns: {}'.format(split_line))
                parm_name = split_line[0]
                if parm_name == 'OUTFILE': # special case because there are multiple occurrences, not consecutive
                        n_outfile_lines += 1
                        parm_name = 'OUTFILE_' + str(n_outfile_lines)
                elif parm_name == 'OUTVAR': # special case because multiple OUTVAR lines follow each OUTFILE line
                        parm_name = 'OUTVAR_' + str(n_outfile_lines)
                try:
                    if global_parms[parm_name]: # if we've already read one or more entries of this parm_name
#                        print('parm {} exists already'.format(parm_name))
                        global_parms[parm_name].append(split_line[1:])
                except:
                    global_parms[parm_name] = []
                    global_parms[parm_name].append(split_line[1:])
#                    print('global_parms[{}]: {}'.format(parm_name,global_parms[parm_name]))
                # We need to create a placeholder in this position for INIT_STATE if it doesn't exist in the initial
                # global parameters file, to be used for all iterations after the first year
                if parm_name == 'OUTPUT_FORCE': # OUTPUT_FORCE should always immediately precede INIT_STATE in the global file
                    global_parms['INIT_STATE'] = []
    return global_parms

def update_global_parms(global_parms, temp_vpf, temp_snb, num_snow_bands, start_year, start_month, start_day, end_year, end_month, end_day, init_state_file, state_save_year, state_save_month, state_save_day):
    """ Updates the global_parms dict at the beginning of an annual iteration """
    # Set start and end dates for the upcoming VIC run (should be one year long, except for the spin-up run)
    global_parms['STARTYEAR'] = [[str(start_year)]]
    global_parms['STARTMONTH'] = [[str(start_month)]]
    global_parms['STARTDAY'] = [[str(start_day)]]
    global_parms['ENDYEAR'] = [[str(end_year)]]
    global_parms['ENDMONTH'] = [[str(end_month)]]
    global_parms['ENDDAY'] = [[str(end_day)]]
    # Set new output state file parameters for upcoming VIC run
    global_parms['STATEYEAR'] = state_save_year
    global_parms['STATEMONTH'] = state_save_month
    global_parms['STATEDAY'] = state_save_day
    global_parms['VEGPARAM'] = [[temp_vpf]]
    global_parms['SNOW_BAND'] = [[num_snow_bands, temp_snb]]
    # All VIC iterations except for the initial spin-up period have to load a saved state from the previous
    if init_state_file:
        global_parms['INIT_STATE'] = [[init_state_file]]
        
def write_global_parms_file(global_parms, temp_gpf):
    """ Reads existing global_parms dict and writes out a new temporary VIC Global Parameter File for feeding into VIC to the temporary global parameters file temp_gpf """
    with open(temp_gpf, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        for parm_name, parm_value in global_parms.items():
            num_parm_lines = len(global_parms[parm])
            if parm == 'INIT_STATE' and len(global_parms['INIT_STATE']) == 0:
                pass
            elif parm[0:8] == 'OUTFILE_':
                line = []
                line.append('OUTFILE')
                for value in parm_value[0]:
                    line.append(value)
                writer.writerow(line)
            elif parm[0:7] == 'OUTVAR_':
                for line_num in range(num_parm_lines):
                    line = []
                    line.append('OUTVAR')
                    for value in parm_value[line_num]:
                        line.append(value)
                        writer.writerow(line)
            elif num_parm_lines == 1:
                line = []
                line.append(parm)
                for value in parm_value[0]:
                    line.append(value)
                writer.writerow(line)
            elif num_parm_lines > 1:
                for line_num in range(num_parm_lines):
                    line = []
                    line.append(parm)
                    for value in parm_value[line_num]:
                        line.append(value)
                    writer.writerow(line)


def get_veg_parms(veg_parm_file):
    """Reads in a Vegetation Parameter File and parses out VIC grid cell IDs,
       as well as an ordered dict of all vegetation parameters,
       grouped by elevation band index
    """
    vp = VegParams(veg_parm_file)
    return vp, vp.cell_ids

def init_residual_area_fracs(cell_ids, veg_parms, snb_parms):
    """ Reads the initial snow band area fractions and glacier vegetation (HRU) tile area fractions and calculates the initial residual area fractions """
    residual_area_fracs = {}
    for cell in cell_ids:
        print('cell: {}'.format(cell))
        residual_area_fracs[cell] = {}
        for band in veg_parms[cell]:
            print('band: {}'.format(band))
            glacier_exists = False
            for line_idx, line in enumerate(veg_parms[cell][band]):
                print('line: {}'.format(line))
                if line[0] == GLACIER_ID:
                    glacier_exists = True
                    residual_area_fracs[cell][band] = snb_parms[cell][0][int(band)] - float(veg_parms[cell][band][line_idx][1])
                    print('Glacier exists in this band. residual_area_fracs[{}][{}] = {} - {} = {}'.format(cell, band, snb_parms[cell][0][int(band)], float(veg_parms[cell][band][line_idx][1]), residual_area_fracs[cell][band]))
                    if residual_area_fracs[cell][band] < 0:
                        print('init_residual_area_fracs(): Error: Calculated a negative residual area fraction for cell {}, band {}. The sum of vegetation tile fraction areas for a given band in the Vegetation Parameter File must be equal to the area fraction for that band in the Snow Band File. Exiting.\n'.format(cell, band))
                        sys.exit(0)
                    break
            if not glacier_exists:
                residual_area_fracs[cell][band] = snb_parms[cell][0][int(band)]
                print('No glacier in this band.  residual_area_fracs[{}][{}] = {}'.format(cell, band, residual_area_fracs[cell][band]))
    return residual_area_fracs

def get_initial_residual_area_fracs():
    """ Reads initial vegetation tile area fractions for all non-glacier
        vegetation types (HRUs)
    """
    pass

def update_band_areas(cell_ids, band_map, num_snow_bands, band_size, pixel_to_cell_map,
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

def update_veg_parms(cell_ids, veg_parms, area_frac_bands, area_frac_glacier, residual_area_fracs):
    """ Updates vegetation parameters for all VIC grid cells by applying calculated changes in glacier area fractions across all elevation bands """
    for cell in cell_ids:
        print('cell: {}'.format(cell))
        for band in veg_parms[cell]:
            if area_frac_bands[cell][int(band)] > 0: # If we (still) have anything in this elevation band
                new_residual_area_frac = area_frac_bands[cell][int(band)] - area_frac_glacier[cell][int(band)]
                delta_residual = residual_area_fracs[cell][band] - new_residual_area_frac
                veg_types = [] # identify and temporarily store vegetation types (HRUs) currently existing within this band
                bare_soil_exists = False # temporary boolean to identify if bare soil vegetation type (HRU) currently exists within this band
                glacier_exists = False # temporary boolean to identify if glacier vegetation type (HRU) currently exists within this band 
                sum_previous_band_area_fracs = 0 # sum of band vegetation tile (HRU) area fractions in the last iteration
                if delta_residual < 0: # glacier portion of this band has SHRUNK; update its area fraction and increase the bare soil component accordingly
                    print('\nGlacier portion of band {} has SHRUNK (delta_residual = {})'.format(band, delta_residual))
                    for line_idx, line in enumerate(veg_parms[cell][band]): 
                        print('Current line in band {}:  {}'.format(band, line))
                        veg_types.append(int(line[0]))
                        print('Existing vegetation type identified in this elevation band: {}'.format(veg_types[-1]))
                        if veg_types[-1] == int(BARE_SOIL_ID):
                            bare_soil_exists = True
                            veg_parms[cell][band][line_idx][1] = str(float(veg_parms[cell][band][line_idx][1]) + abs(delta_residual))
                            print('Bare soil already exists in this band.  New area fraction: veg_parms[{}][{}][{}][1] = {}'.format(cell, band, line_idx, veg_parms[cell][band][line_idx][1]))
                        elif veg_types[-1] == int(GLACIER_ID):
                            veg_parms[cell][band][line_idx][1] = str(area_frac_glacier[cell][int(band)])
                            print('Updated glacier area fraction: veg_parms[{}][{}][{}][1] = {}'.format(cell, band, line_idx, veg_parms[cell][band][line_idx][1]))
                    if not bare_soil_exists: # insert new line for bare soil vegetation (HRU) type in the correct numerical position within this band
                        #new_line = BARE_SOIL_ID + ' ' + str(delta_residual) + ' ' + ''.join(str(x) for x in bare_soil_root_parms) + ' ' + str(band)
                        new_line = []
                        new_line.append(BARE_SOIL_ID)
                        new_line.append(str(delta_residual))
                        for num in bare_soil_root_parms:
                            new_line.append(str(num))
                        new_line.append(str(band))                        
                        bisect.insort_left(veg_types, int(BARE_SOIL_ID))
                        position = veg_types.index(int(BARE_SOIL_ID))
                        veg_parms[cell][band].insert(position, new_line)
                        print('Bare soil did NOT already exist in this band. Inserted line: veg_parms[{}][{}][{}] = {}'.format(cell, band, position, veg_parms[cell][band][position]))
                elif delta_residual > 0: # glacier portion of this band has GROWN; decrease fraction of non-glacier HRUs by a total of delta_residual
                    print('\nGlacier portion of band {} has GROWN (delta_residual = {})'.format(band, delta_residual))
                    for line_idx, line in enumerate(veg_parms[cell][band]): 
                        sum_previous_band_area_fracs += float(line[1])
                    print('sum_previous_band_area_fracs = {}'.format(sum_previous_band_area_fracs))
                    for line_idx, line in enumerate(veg_parms[cell][band]):
                        print('Current line in band {}:  {}'.format(band, line))
                        veg_types.append(int(line[0])) # keep a record of vegetation types (HRUs) currently existing within this band
                        print('Existing vegetation type identified in this elevation band: {}'.format(veg_types[-1]))
                        if veg_types[-1] == int(GLACIER_ID): # set glacier tile area fraction
                            glacier_exists = True
                            veg_parms[cell][band][line_idx][1] = str(area_frac_glacier[cell][int(band)])
                            print('Glacier already exists in this band.  New glacier area fraction: veg_parms[{}][{}][{}][1] = {}'.format(cell, band, line_idx, veg_parms[cell][band][line_idx][1]))
                        else:
                            # Get the change in area fraction for this vegetation type (HRU) based on its previous share of the residuals in this elevation band
                            print('Existing area fraction for this vegetation type: {}'.format(veg_parms[cell][band][line_idx][1]))
                            delta_veg_area_frac = delta_residual * (float(veg_parms[cell][band][line_idx][1]) / sum_previous_band_area_fracs)
                            print('Calculated delta_veg_area_frac = {}'.format(delta_veg_area_frac))
                            new_veg_area_frac = float(veg_parms[cell][band][line_idx][1]) - delta_veg_area_frac
                            if new_veg_area_frac >= 0:
                                veg_parms[cell][band][line_idx][1] = str(new_veg_area_frac)
                                print('Reduced area fraction of vegetation type {} by delta_veg_area_frac.  New area fraction: veg_parms[{}][{}][{}][1] = {}'.format(veg_types[-1], cell, band, line_idx, veg_parms[cell][band][line_idx][1]))
                            else: #the existing vegetation tile (HRU) was overtaken by glacier
                                # delete the existing vegetation tile (HRU) that was overtaken by glacier?
                                #print('Area fraction of vegetation type {} in cell {}, band {} was reduced to {}. Removing tile (HRU) from vegetation parameters'.format(veg_types[-1], cell, band, new_veg_area_frac)
                                #del veg_parms[cell][band][line_idx]
                                print('update_veg_parms(): Error: Calculated a negative area fraction ({}) for vegetation type {} in cell {}, band {}. Exiting.\n'.format(new_veg_area_frac, veg_types[-1], cell, band))
                                sys.exit(0)
                    if not glacier_exists: # insert new line for glacier vegetation (HRU) type in the correct numerical position within this band
                        #new_line = GLACIER_ID + ' ' + str(area_frac_glacier[cell][int(band)]) + ' ' + ''.join(str(x) for x in glacier_root_parms) + ' ' + str(band)
                        new_line = []
                        new_line.append(GLACIER_ID)
                        new_line.append(str(area_frac_glacier[cell][int(band)]))
                        for num in glacier_root_parms:
                            new_line.append(str(num))
                        new_line.append(str(band))
                        bisect.insort_left(veg_types, int(GLACIER_ID))
                        position = veg_types.index(int(GLACIER_ID))
                        veg_parms[cell][band].insert(position, new_line)
                        print('Glacier did NOT already exist in this band. Inserted line with new glacier area fraction: veg_parms[{}][{}][{}] = {}'.format(cell, band, position, veg_parms[cell][band][position]))
            else: # We have nothing left in this elevation band.  Set all fractional areas for all vegetation tiles (HRUs) in this band to zero
                print('\nNothing left in elevation band {}.  Setting all fractional areas to zero.'.format(band))
                for line_idx, line in enumerate(veg_parms[cell][band]):
                    veg_parms[cell][band][line_idx][1] = 0
            # Set residual area fractions to the new calculated values for the next iteration
            residual_area_fracs[cell][band] = new_residual_area_frac

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

def create_band_map(snb_parms, band_size):
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
    for cell in snb_parms:
        #print('cell: {}'.format(cell)
        snb_parms[cell][0] = area_frac_bands[cell]
        print(' '.join(map(str, snb_parms[cell][0])))

def write_snb_parms_file(temp_snb, snb_parms, area_frac_bands):
    """ Writes current (updated) snow band parameters to a new temporary
        Snow Band File for feeding back into VIC
    """
    with open(temp_snb, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        for cell in snb_parms:
            line = [ cell +
                     area_frac_bands[cell] +
                     snb_parms[cell][1] + # median elevations
                     snb_parms[cell][2] # Pfactor values
                   ]
            writer.writerow(line)

# Main program.  
if __name__ == '__main__':
    main()

def main():
    print('\n\nVIC + RGM ... together at last!')

    # Parse command line parameters
    vic_global_file, rgm_params_file, surf_dem_in_file, bed_dem_file, pixel_cell_map_file, \
    init_glacier_mask_file, output_trace_files, glacier_root_parms_file, band_size = parse_input_parms()

    # Get all initial VIC global parameters from the global parameter file
    global_parms = get_global_parms(vic_global_file)

    # Get entire time range of coupled VIC-RGM run from the initial VIC global file
    start_year = int(global_parms['STARTYEAR'][0][0])
    start_month = int(global_parms['STARTMONTH'][0][0])
    start_day = int(global_parms['STARTDAY'][0][0])
    start_date = date(start_year, start_month, start_day)
    end_year = int(global_parms['ENDYEAR'][0][0])
    end_month = int(global_parms['ENDMONTH'][0][0])
    end_day = int(global_parms['ENDDAY'][0][0])
    end_date = date(end_year, end_month, end_day)
    # Set the initial year for the coupled VIC-RGM simulation
    year = start_date.year
    # Get the date that glacier accumulation is to start (at the end of VIC "spin-up")
    glacier_accum_start_year = int(global_parms['GLACIER_ACCUM_START_YEAR'][0][0])
    glacier_accum_start_month = int(global_parms['GLACIER_ACCUM_START_MONTH'][0][0])
    glacier_accum_start_day = int(global_parms['GLACIER_ACCUM_START_DAY'][0][0])
    glacier_accum_start_date = date(glacier_accum_start_year, glacier_accum_start_month, glacier_accum_start_day)

    # Initial VIC output state filename prefix is determined by STATENAME in the global file
    state_filename_prefix = global_parms['STATENAME'][0][0]
    # Numeric code indicating a glacier vegetation tile (HRU)
    GLACIER_ID = global_parms['GLACIER_ID'][0]

    # Get VIC vegetation parameters and grid cell IDs from initial Vegetation Parameter File
    veg_parm_file = global_parms['VEGPARAM'][0][0]
    veg_parms, cell_ids = get_veg_parms(veg_parm_file)

    # Get VIC snow/elevation band parameters from initial Snow Band File
    num_snow_bands = int(global_parms['SNOW_BAND'][0][0])
    snb_file = global_parms['SNOW_BAND'][0][1]
    snb_parms = get_snb_parms(snb_file, num_snow_bands)

    # Get list of elevation bands for each VIC grid cell
    band_map = create_band_map(snb_parms, band_size)

    # The RGM will always output a DEM file of the same name (if running RGM for a single year at a time)
    rgm_surf_dem_out_file = temp_files_path + 's_out_00001.grd'

    # Calculate the initial residual (i.e. non-glacier) area fractions for all bands in all cells 
    #residual_area_fracs = init_residual_area_fracs(cell_ids, veg_parms, snb_parms)
    # NOTE: this should probably go move to after the initial update_band_areas() and before update_veg_parms()

    # Open and read VIC-grid-to-RGM-pixel mapping file
    # pixel_to_cell_map is a list of dimensions num_rows_dem x num_cols_dem, each element containing a VIC grid cell ID
    pixel_to_cell_map, num_rows_dem, num_cols_dem, cell_areas = get_rgm_pixel_mapping(pixel_cell_map_file)

    # Get DEM xmin, xmax, ymin, ymax metadata of Bed DEM and check file header validity     
    dem_xmin, dem_xmax, dem_ymin, dem_ymax, num_rows, num_cols = read_gsa_headers(bed_dem_file)
    # Verify number of columns & rows agree with what's stated in the pixel_to_cell_map_file
    assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch of stated dimension(s) \
        between Bed DEM in {} (num rows: {}, num columns: {}) and the VIC-grid-to-RGM-pixel map in {} \
        (num rows: {}, num columns: {}). Exiting.\n'.format(bed_dem_file, num_rows, num_cols, pixel_cell_map_file, num_rows_dem, num_cols_dem))
    # Read in the provided Bed Digital Elevation Map (BDEM) file to 2D bed_dem array
    bed_dem = np.loadtxt(bed_dem_file, skiprows=5)

    # Check header validity of Surface DEM file
    _, _, _, _, num_rows, num_cols = read_gsa_headers(surf_dem_in_file)
    # Verify number of columns & rows agree with what's stated in the pixel_to_cell_map_file
    assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch of stated dimension(s) \
        between Surface DEM in {} (num rows: {}, num columns: {}) and the VIC-grid-to-RGM-pixel map in {} \
        (num rows: {}, num columns: {}). Exiting.\n'.format(surf_dem_in_file, num_rows, num_cols, pixel_cell_map_file, num_rows_dem, num_cols_dem))
    # Read in the provided Surface Digital Elevation Map (SDEM) file to 2D surf_dem array
    surf_dem_initial = np.loadtxt(surf_dem_in_file, skiprows=5)

    # Check header validity of Glacier Mask file
    _, _, _, _, num_rows, num_cols = read_gsa_headers(glacier_mask_file)
    # Verify number of columns & rows agree with what's stated in the pixel_to_cell_map_file
    assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch of stated dimension(s) \
        between Glacier Mask in {} (num rows: {}, num columns: {}) and the VIC-grid-to-RGM-pixel map in {} \
        (num rows: {}, num columns: {}). Exiting.\n'.format(glacier_mask_file, num_rows, num_cols, pixel_cell_map_file, num_rows_dem, num_cols_dem))
    # Read in the provided initial glacier mask file to 2D glacier_mask array 
    glacier_mask = np.loadtxt(glacier_mask_file, skiprows=5)

    # Apply the initial glacier mask and modify the band and glacier area fractions accordingly
    area_frac_bands, area_frac_glacier = update_band_areas(cell_ids, band_map, num_snow_bands, band_size, pixel_to_cell_map,
                      surf_dem_initial, num_rows_dem, num_cols_dem)

    # Calculate the initial residual (i.e. non-glacier) area fractions for all bands in all cells
    residual_area_fracs = init_residual_area_fracs(cell_ids, veg_parms, snb_parms)

    # Update the vegetation parameters vis-a-vis the application of the initial glacier mask, and write to new temporary file temp_vpf
    update_veg_parms(cell_ids, veg_parms, area_frac_bands, area_frac_glacier, residual_area_fracs)
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
        #    start_year, start_month, start_day, \
        #    end_year, end_month, end_day, \
        #    init_state_file, state_save_year, state_save_month, state_save_day)
        if year == start_date.year:
            # set global parameters for VIC "spin-up", from start_date.year to the end of the first glacier accumulation
            init_state_file = None
            temp_end_date = glacier_accum_start_date - timedelta(days=1)
            update_global_parms(global_parms, temp_vpf, temp_snb, num_snow_bands, \
                start_date.year, start_date.month, start_date.day, \
                temp_end_date.year, temp_end_date.month, temp_end_date.day, \
                init_state_file, temp_end_date.year, temp_end_date.month, temp_end_date.day)
            # Fast-forward year to what it will be when VIC finishes its spin-up
            year = temp_end_date.year
        else:
            temp_start_date = date(year, glacier_accum_start_date.month, glacier_accum_start_date.day)
            temp_end_date = date(year+1, glacier_accum_start_date.month, glacier_accum_start_date.day) - timedelta(days=1)
            # set/create INIT_STATE parm as the last written VIC state_file (parm does not exist in the first read-in of global_parms)
            temp_state_date = temp_end_date - timedelta(days=1)
            init_state_file = state_filename_prefix + "_" + str(year-1) + str(temp_state_date.month) + str(temp_state_date.day)
            update_global_parms(global_parms, temp_vpf, temp_snb, num_snow_bands, \
                temp_start_date.year, temp_start_date.month, temp_start_date.day, \
                temp_end_date.year, temp_end_date.month, temp_end_date.day, \
                init_state_file, temp_end_date.year, temp_end_date.month, temp_end_date.day)

        write_global_parms_file(global_parms, temp_gpf)
        print('invoking VIC with global parameter file {}'.format(temp_gpf))

        # 2. Run VIC for a year.  This will save VIC model state at the end of the year, along with a Glacier Mass Balance (GMB) polynomial for each cell
        subprocess.check_call([vic_full_path, "-g", temp_gpf], shell=False, stderr=subprocess.STDOUT)

        # 3. Open VIC NetCDF state file and get the most recent GMB polynomial for each grid cell being modeled
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
        area_frac_bands, area_frac_glacier = update_band_areas(cell_ids, band_map, num_snow_bands, band_size, pixel_to_cell_map, rgm_surf_dem_out, num_rows_dem, num_cols_dem)

        # 9. Update vegetation parameters and write to new temporary file temp_vpf
        update_veg_parms(cell_ids, veg_parms, area_frac_bands, area_frac_glacier, residual_area_fracs)
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
