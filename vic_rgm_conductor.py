#!/usr/bin/env python

""" This script orchestrates a coupled Variable Infiltration Capacity (VIC) and Regional Glacier Model (RGM) run. """

import sys
import os
import argparse
import subprocess
import numpy as np
import h5py
import csv
from collections import OrderedDict
#import collections

vic_full_path = '/home/mfischer/code/vic/vicNl'
rgm_full_path = '/home/mfischer/code/rgm/rgm'
input_files_path = '/home/mfischer/vic_dev/input/peyto/'
#rgm_params_file =  input_files_path + 'global_params_VIC.txt'
#initial_vic_global_file = '/home/mfischer/vic_dev/input/place/glb_base_PLACE_19601995_VIC4.1.2_outNETCDF_initial.txt' 
#vic_global_file = '/home/mfischer/vic_dev/input/place/glb_base_PLACE_19601995_VIC4.1.2_outNETCDF.txt' 
#vpf_full_path = '/home/mfischer/vic_dev/input/place/vpf_place_100m.txt' 
#rgm_output_path = '/home/mfischer/vic_dev/out/testing/rgm_output/'
temp_files_path = '/home/mfischer/vic_dev/out/testing/temp_out_files/'
# NOTE: Setting a default elevation band size of 100 m (should this be a command line parameter?)
band_size = 100

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s]n' % message)
        self.print_help()
        sys.exit(2)

# Get all global parameters 
parser = MyParser()
parser.add_argument('--g', action="store", dest="vic_global_file", type=str, help = 'file name and path of the VIC global parameters file')
parser.add_argument('--rgm-params', action="store", dest="rgm_params_file", type=str, help = 'file name and path of the Regional Glacier Model (RGM) parameters file')
parser.add_argument('--sdem', action="store", dest="surf_dem_file", type=str, help = 'file name and path of the initial Surface Digital Elevation Model (SDEM) file')
parser.add_argument('--bdem', action="store", dest="bed_dem_file", type=str, help = 'file name and path of the Bed Digital Elevation Model (BDEM) file')
parser.add_argument('--pixel-map', action="store", dest="pixel_cell_map_file", type=str, help = 'file name and path of the RGM Pixel to VIC Grid Cell mapping file')
parser.add_argument('--trace-files', action="store_true", default=False, dest="trace_files", help = 'write out persistent ASCII DEM and mass balance grid files on each time step for offline inspection')

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
options = parser.parse_args()
vic_global_file = options.vic_global_file
rgm_params_file = options.rgm_params_file
rgm_surf_dem_in_file = options.surf_dem_file
bed_dem_file = options.bed_dem_file
pixel_cell_map_file = options.pixel_cell_map_file
output_trace_files = options.trace_files


# To have nested ordered defaultdicts
class OrderedDefaultdict(OrderedDict):
	# from: http://stackoverflow.com/questions/4126348/how-do-i-rewrite-this-function-to-implement-ordereddict/4127426#4127426
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
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
        return self.__class__, args, None, None, self.iteritems()

def get_global_parms(global_parm_file):
	""" Parses the initial VIC global parameters file created by the user with the settings for the entire VIC-RGM run """
	global_parms = OrderedDefaultdict()
	n_outfile_lines = 0
	with open(global_parm_file, 'r') as f:
		for line in iter(f):
			#print 'line: {}'.format(line)
			if not line.isspace() and line[0] is not '#':
				split_line = line.split()
				#print 'columns: {}'.format(split_line)
				parm_name = split_line[0]
				if parm_name == 'OUTFILE': # special case because there are multiple occurrences, not consecutive
						n_outfile_lines += 1
						parm_name = 'OUTFILE_' + str(n_outfile_lines)
				elif parm_name == 'OUTVAR': # special case because multiple OUTVAR lines follow each OUTFILE line
						parm_name = 'OUTVAR_' + str(n_outfile_lines)
				try:
					if global_parms[parm_name]: # if we've already read one or more entries of this parm_name
#						print 'parm {} exists already'.format(parm_name)
						global_parms[parm_name].append(split_line[1:])
				except:
					global_parms[parm_name] = []
					global_parms[parm_name].append(split_line[1:])
#					print 'global_parms[{}]: {}'.format(parm_name,global_parms[parm_name])
				# We need to create a placeholder in this position for INIT_STATE if it doesn't exist in the initial
				# global parameters file, to be used for all iterations after the first year
				if parm_name == 'OUTPUT_FORCE': # OUTPUT_FORCE should always immediately precede INIT_STATE in the global file
					global_parms['INIT_STATE'] = []
	return global_parms

def update_global_parms(init_state):
	# Important parms: STARTYEAR, ENDYEAR, VEGPARAM, SNOW_BAND (and GLACIER_ACCUM_START_YEAR, GLACIER_ACCUM_INTERVAL?)
	global_parms['STARTYEAR'] = str(year)
	global_parms['ENDYEAR'] = str(year)
	# All iterations after the first / wind-up period have modified state, vegetation parms, and snow band parms
	if init_state:
		# set/create INIT_STATE parm with most current state_file (does not exist in the first read-in of global_parms)
		init_state_file = state_filename_prefix + "_" + str(year - 1) + str(global_parms['STATEMONTH'][0]) + str(global_parms['STATEDAY'][0])
		global_parms['INIT_STATE'].append(init_state_file)
		# New output state filename for next VIC year run
		global_parms['STATENAME'] = [state_file]
		global_parms['VEGPARAM'] = temp_vpf
		global_parms['SNOW_BAND'] = temp_snb

def write_global_parms_file():
	""" Reads existing global_parms dict and writes out a new temporary VIC Global Parameter File for feeding into VIC """
	temp_gpf = temp_files_path + 'gpf_temp_' + str(year) + '.txt'
	with open(temp_gpf, 'w') as f:
		writer = csv.writer(f, delimiter=' ')
		for parm in global_parms:
			num_parm_lines = len(global_parms[parm])
			if parm == 'INIT_STATE' and len(global_parms['INIT_STATE']) == 0:
				pass
			elif parm[0:8] == 'OUTFILE_':
				line = []
				line.append('OUTFILE')
				for value in global_parms[parm][0]:
					line.append(value)
				writer.writerow(line)
			elif parm[0:7] == 'OUTVAR_':
				for line_num in range(0, num_parm_lines):
					line = []
					line.append('OUTVAR')
					for value in global_parms[parm][line_num]:
						line.append(value)
						writer.writerow(line)
			elif num_parm_lines == 1:
				line = []
				line.append(parm)
				for value in global_parms[parm][0]:
					line.append(value)
				writer.writerow(line)
			elif num_parm_lines > 1:
				for line_num in range(0, num_parm_lines):
					line = []
					line.append(parm)
					for value in global_parms[parm][line_num]:
						line.append(value)
					writer.writerow(line)
	return temp_gpf

def get_veg_parms(veg_parm_file):
	""" Reads in a Vegetation Parameter File and parses out VIC grid cell IDs, as well as an ordered nested dict of all vegetation parameters,
	grouped by elevation band index """
	cell_ids = []
	num_veg_tiles = {}
	veg_parms = OrderedDefaultdict()
	with open(veg_parm_file, 'r') as f:
		for line in iter(f):
			split_line = line.split()
			num_columns = len(split_line)
			if num_columns == 2:
				cell = split_line[0]
				cell_ids.append(cell)
				num_veg_tiles[cell] = split_line[1]
				veg_parms[cell] = OrderedDefaultdict()
			elif num_columns > 2:
				band_id = split_line[-1]
				try:
					veg_parms[cell][band_id].append(split_line)
				except:
					veg_parms[cell][band_id] = []
					veg_parms[cell][band_id].append(split_line)
	return cell_ids, num_veg_tiles, veg_parms

def update_veg_parms():
	""" Updates vegetation parameters for all VIC grid cells by applying calculated changes in glacier area fractions across all elevation bands """
	GLACIER_ID = global_parms['GLACIER_ID'][0]
	for cell in cell_ids:
#		print 'cell: {}'.format(cell)
		for band in veg_parms[cell]:
			num_tiles_in_band = len(veg_parms[cell][band])
#			print 'band: {}, num_tiles_in_band: {}'.format(band, num_tiles_in_band)
			if area_frac_glacier[cell][int(band)] > 0: # there exists a glacier tile in this band
#				print 'glacier exists in band {}'.format(band)
				# the remaining area fraction to be distributed among non-glacier tiles:
				area_frac_to_distribute = 1 - area_frac_glacier[cell][int(band)]
#				print 'area_frac_to_distribute: {}'.format(area_frac_to_distribute)
				# the updated evenly-distributed area fraction values for each non-glacier tile: 
				non_glacier_single_tile_area_frac = area_frac_to_distribute / (num_tiles_in_band - 1)
#				print 'non_glacier_single_tile_area_frac: {}'.format(non_glacier_single_tile_area_frac)
				for line_idx, line in enumerate(veg_parms[cell][band]): # go through all tile lines in this band
#					print 'line in band {}:  {}'.format(band, line)
					if line[0] == GLACIER_ID: 
						veg_parms[cell][band][line_idx][1] = str(area_frac_glacier[cell][int(band)])
					else:
						veg_parms[cell][band][line_idx][1] = str(non_glacier_single_tile_area_frac)

def write_veg_parms_file():
	""" Writes current (updated) vegetation parameters to a new temporary Vegetation Parameters File for feeding back into VIC """
	temp_vpf = temp_files_path + 'vpf_temp_' + str(year) + '.txt'
	with open(temp_vpf, 'w') as f:
		writer = csv.writer(f, delimiter=' ')
		for cell in veg_parms:
			writer.writerow([cell, num_veg_tiles[cell]])
			for band in veg_parms[cell]:
				for line in veg_parms[cell][band]:
					writer.writerow(line)
	return temp_vpf

def get_snb_parms(snb_file, num_snow_bands):
	""" Reads in a Snow Band File and outputs an ordered dict:
	{'cell_id_0' : [area_frac_band_0,...,area_frac_band_N],[median_elev_band_0,...,median_elev_band_N],[Pfactor_band_0,...,Pfactor_band_N]], 'cell_id_1' : ..."""
	snb_parms = OrderedDict()
	with open(snb_file, 'r') as f:
		for line in iter(f):
			#print 'snb file line: {}'.format(line)
			split_line = line.split()
			num_columns = len(split_line)
			if num_columns != 3*num_snow_bands + 1:
				print 'get_snb_parms(): Number of columns ({}) in snow band file {} is incorrect for the given number of snow bands ({}) given in the global parameter file (should be 3 * num_snow_bands + 1). Exiting.\n'.format(num_columns, snb_file, num_snow_bands)
				#sys.exit(0)
			snb_parms[split_line[0]] = [[float(x) for x in split_line[1 : num_snow_bands+1]],[float(x) for x in split_line[num_snow_bands+1 : 2*num_snow_bands+1]],[float(x) for x in split_line[2*num_snow_bands+1 : 3*num_snow_bands+1]]]
	return snb_parms

def update_snb_parms():
	for cell in snb_parms:
		snb_parms[cell][0] = area_frac_bands[cell]

def write_snb_parms_file():
	""" Writes current (updated) snow band parameters to a new temporary Snow Band File for feeding back into VIC """
	temp_snb = temp_files_path + 'snb_temp_' + str(year) + '.txt'
	with open(temp_snb, 'w') as f:
		writer = csv.writer(f, delimiter=' ')
		for cell in snb_parms:
			line = []
			line.append(cell)
			for area_frac in area_frac_band[cell]:
				line.append(area_frac)
			for band_frac in snb_parms[cell][1]:
				line.append(band_frac) # append existing median elevations
			for pfactor in snb_parms[cell][2]:
				line.append(pfactor) # append existing Pfactor values
			writer.writerow(line)
	return temp_snb

def get_bands(snb_parms, band_size):
	""" Takes a dict of Snow Band parameters and identifies and creates a list of elevation bands for each grid cell of band_size width in meters"""
	band_map = {}
	for cell in cell_ids:
		band_map[cell] = []
		for band in snb_parms[cell][1]:
			if band > 0: # ignore bands with median elevation 0
				band_map[cell].append(int(band - band % band_size))
	return band_map

def get_rgm_pixel_mapping(pixel_map_file):
	""" Parses the RGM pixel to VIC grid cell mapping file and initialises a 2D  
	   	grid of dimensions num_rows_rpg x num_cols_rpg (matching the RGM pixel grid),
	   	each element containing the VIC cell ID associated with that RGM pixel"""
	pixel_grid = []
	with open(pixel_map_file, 'r') as f:
		num_cols_rpg = 0
		num_rows_rpg = 0
		for line in iter(f):
			#print 'line: {}'.format(line)
			split_line = line.split()
			if split_line[0] == 'NCOLS':
				num_cols_rpg = int(split_line[1])
				#print 'num_cols_rpg: {}'.format(num_cols_rpg)
				if num_cols_rpg and num_rows_rpg:
					pixel_grid = [[0 for x in range(0, num_cols_rpg)] for x in range(0, num_rows_rpg)]
			elif split_line[0] == 'NROWS':
				num_rows_rpg = int(split_line[1])
				#print 'num_rows_rpg: {}'.format(num_rows_rpg)
				if num_cols_rpg and num_rows_rpg:
					pixel_grid = [[0 for x in range(0, num_cols_rpg)] for x in range(0, num_rows_rpg)]
			elif split_line[0][0] == '"': # column header row / comments
				#print 'comment line: {}'.format(split_line)
				pass
			else:
				# NOTE: we might want Markus to recreate this mapping file with zero-based indexing
				row_num = int(split_line[1])
				col_num = int(split_line[2])
				#print 'populating row: {}, col: {}'.format(row_num, col_num)
				# NOTE: changed to only putting VIC cell ID into the row,col position (band and elev are not static)
				pixel_grid[row_num-1][col_num-1] = split_line[5]
			#print 'columns: {}'.format(columns)
	return pixel_grid, num_rows_rpg, num_cols_rpg

def read_asc_dem_file(dem_file):
	""" Opens a DEM file, parses RGM grid extents parameters, and reads in the DEM """
	# NOTE: this may have to change to conform to ARC/ASCII grid format instead
	dem_grid = []
	with open(dem_file, 'r') as f:
		num_cols_dem = 0
		num_rows_dem = 0
		for line in iter(f):
			split_line = line.split()
			if split_line[0] == 'NCOLS':
				num_cols_dem = int(split_line[1])
				#print 'num_cols_dem: {}'.format(num_cols_dem)
				# we assume here that get_rgm_pixel_mapping() has already been called
				if num_cols_dem != num_cols_rpg:
					print 'read_asc_dem_file({}): NCOLS ({}) disagrees with NCOLS ({}) stated in RGM-VIC mapping file. Exiting.\n'.format(dem_file, num_cols_dem, num_cols_rpg)
					sys.exit(0)
			elif split_line[0] == 'NROWS':
				num_rows_dem = int(split_line[1])
				#print 'num_rows_dem: {}'.format(num_rows_bdem)
				# we assume here that get_rgm_pixel_mapping() has already been called
				if num_rows_dem != num_rows_rpg:
					print 'read_asc_dem_file({}): NROWS ({}) disagrees with NROWS ({}) stated in RGM-VIC mapping file. Exiting.\n'.format(dem_file, num_rows_dem, num_rows_rpg)
					sys.exit(0)
			elif split_line[0] == 'XLLCORNER':
				xmin = float(split_line[1])
				#print 'xmin: {}'.format(xmin)
			elif split_line[0] == 'YLLCORNER':
				ymin = float(split_line[1])
				#print 'ymin: {}'.format(ymin)
			elif split_line[0] == 'CELLSIZE':
				dem_pixel_size = int(split_line[1])
				# Calculate the pixel area in square meters
				dem_pixel_area = dem_pixel_size * dem_pixel_size
				if xmin and ymin and dem_pixel_size:
					xmax = (num_cols_dem - 1) * dem_pixel_size + xmin
					ymax = (num_rows_dem - 1) * dem_pixel_size + ymin
					#print 'xmax: {}'.format(xmax)
					#print 'ymax: {}'.format(ymax)
			elif split_line[0] == 'NODATA_value': # should we do anything with this?
				dem_no_data_value = float(split_line[1])
			else: # read the next row of DEM data in
				split_line = [float(x) for x in line.split()] # bad that we're doing the split() operation twice, but this fcn is probably only temporary anyway 
				dem_grid.append(split_line)
	return dem_grid, dem_pixel_area, xmin, xmax, ymin, ymax

def read_gsa_dem_file(dem_file):
	""" Opens an ASCII grid (.grd) DEM file, parses grid extents parameters, and reads in the DEM """
	dem_grid = []
	with open(dem_file, 'r') as f:
		num_cols_dem = 0
		num_rows_dem = 0
		for line_num, line in enumerate(f):
			if line_num == 0:
				split_line = line.split()
				if split_line[0] != 'DSAA':
					print 'read_gsa_dem_file({}): DSAA header on first line of DEM file was not found or is malformed.  DEM file does not conform to ASCII grid format.  Exiting. \n'.format(dem_file)
					#sys.exit(0)
			elif line_num == 1:
				split_line = [int(x) for x in line.split()]
				num_cols_dem = split_line[0]
				num_rows_dem = split_line[1]
				if (num_cols_dem != num_cols_rpg) or (num_rows_dem != num_rows_rpg):
					print 'read_gsa_dem_file({}): Disagreement in row/column dimensions between DEM file (NROWS={}, NCOLS={}) and RGM-VIC mapping file (NROWS={}, NCOLS={}). Exiting.\n'.format(dem_file, num_rows_dem, num_cols_dem, num_rows_rpg, num_cols_rpg)
					#sys.exit(0)
				#print 'num_rows_dem: {} num_cols_dem: {}'.format(num_rows_dem, num_cols_dem)
			elif line_num == 2:
				split_line = [float(x) for x in line.split()]
				xmin = split_line[0]
				xmax = split_line[1]
			elif line_num == 3:
				split_line = [float(x) for x in line.split()]
				ymin = split_line[0]
				ymax = split_line[1]
				dem_pixel_size_x = (xmax - xmin) / num_cols_dem
				dem_pixel_size_y = (ymax - ymin) / num_rows_dem
				# Ask Markus what to do about disagreement between x & y pixel sizes, and between these and other DEMs in use
#				if dem_pixel_size_x != dem_pixel_size_y: # perhaps these can be close within a number of decimal places?
#					print 'read_gsa_dem_file({}): Calculated DEM x pixel size ({}) does not equal y pixel size ({}). Exiting.\n '.format(dem_file, dem_pixel_size_x, dem_pixel_size_y)
					#sys.exit(0)
				# Calculate the pixel area in square meters
				dem_pixel_area = dem_pixel_size_x * dem_pixel_size_y
			elif line_num == 4:
				pass  # is anything to be done with these zmin, zmax values?
			else: # read the one big row of DEM data into num_rows_bdem & num_cols_bdem list of lists
				# for row_index in range(0, num_rows_dem):
				# 	start_row = row_index * num_cols_dem
				# 	end_row = start_row + num_cols_dem
				# 	row = line[start_row : end_row]
				# 	split_row = [float(x) for x in row.split()]
				# 	print split_row
				# 	dem_grid.append(split_row)
				split_line = [float(x) for x in line.split()] 
				dem_grid.append(split_line)
	return dem_grid, dem_pixel_area, xmin, xmax, ymin, ymax

def get_mass_balance_polynomials(state, state_file):
	""" Extracts the Glacier Mass Balance polynomial for each grid cell from an open VIC state file """
	gmb_info = state['GLAC_MASS_BALANCE_INFO'][0]
	cell_count = len(gmb_info)
	if cell_count != len(cell_ids):
		print 'get_mass_balance_polynomials: The number of VIC cells ({}) read from the state file {} and those read from the vegetation parameter file ({}) disagree. Exiting.\n'.format(cell_count, state_file, len(cell_ids))
		sys.exit(0)
	gmb_polys = {}
	for i in range(0, cell_count):
		cell_id = str(int(gmb_info[i][0]))
		if cell_id not in cell_ids:
			print 'get_mass_balance_polynomials: Cell ID {} was not found in the list of VIC cell IDs read from the vegetation parameters file. Exiting.\n'.format(cell_id)
			sys.exit(0)
		gmb_polys[cell_id] = [gmb_info[i][1], gmb_info[i][2], gmb_info[i][3]]
	return gmb_polys

def mass_balances_to_rgm_grid(gmb_polys, pixel_to_cell_map):
	""" Translate mass balances from grid cell GMB polynomials to 2D RGM pixel grid to use as one of the inputs to RGM """
	mass_balance_grid = [[0 for x in range(0, num_cols_rpg)] for x in range(0, num_rows_rpg)]
	for row in range(0, num_rows_rpg):
		for col in range(0, num_cols_rpg):
			pixel = pixel_to_cell_map[row][col]
			band = pixel[0]
			median_elev = float(pixel[1])
			cell_id = pixel[2]
			# only grabbing pixels that fall within a VIC cell
			if cell_id != 'NA':
				# check that the cell_id agrees with what was read from the veg_parm_file
				if cell_id not in cell_ids:
					print 'mass_balances_to_rgm_grid: Cell ID {} was not found in the list of VIC cell IDs read from the vegetation parameters file. Exiting.\n'.format(cell_id)
					#sys.exit(0)
				mass_balance_grid[row][col] = gmb_polys[cell_id][0] + median_elev * (gmb_polys[cell_id][1] + median_elev * gmb_polys[cell_id][2])
	return mass_balance_grid			
		
def write_grid_to_gsa_file(grid, outfilename):
	""" Writes a 2D grid to ASCII file in the input format expected by the RGM for DEM and mass balance grids """
	zmin = np.min(grid)
	zmax = np.max(grid)
	header_rows = [['DSAA'], [num_cols_rpg, num_rows_rpg], [rgm_xmin, rgm_xmax], [rgm_ymin, rgm_ymax], [zmin, zmax]]
	with open(outfilename, 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter=' ')
		for header in range(0, len(header_rows)):
			writer.writerow(header_rows[header])
		for row in range(0, num_rows_rpg):
			writer.writerow(grid[row])

def update_glacier_mask(sdem, bdem):
	""" Takes output Surface DEM from RGM and uses element-wise differencing with the Bed DEM to form an updated glacier mask """
	diffs = np.array(sdem) - np.array(bed_dem)
	if np.any(diffs < 0):
		print 'update_glacier_mask: Error: subtraction of Bed DEM from output Surface DEM of RGM produced one or more negative values.  Exiting.\n'
		sys.exit(0)
	glacier_mask = np.zeros((num_rows_rpg, num_cols_rpg))
	glacier_mask[diffs > 0] = 1
	glacier_mask[diffs == 0] = 0
	return glacier_mask

def update_band_areas():
	""" Calculates the area fractions of elevation bands within VIC cells, and of glaciers within these bands """
	total_cell_area = num_rows_rpg * num_cols_rpg
	band_areas = {}
	glacier_areas = {}
	area_frac_band = {}
	area_frac_glacier = {}
	for cell in cell_ids:
		band_areas[cell] = []
		glacier_areas[cell] = []
		area_frac_band[cell] = []
		area_frac_glacier[cell] = []
		for band_idx, band in enumerate(band_map[cell]):
			band_areas[cell].append(0)
			glacier_areas[cell].append(0)
			area_frac_band[cell].append(0)
			area_frac_glacier[cell].append(0)
			for row in range(0, num_rows_rpg):
				for col in range(0, num_cols_rpg):
					if pixel_to_cell_map[row][col] == cell: # if this pixel is within the current VIC cell
						if rgm_surf_dem_out[row][col] in range(band, band + band_size-1):
							band_areas[cell][band_idx] += 1
							if glacier_mask[row][col]:
								glacier_areas[cell][band_idx] += 1
			# This will be used to update the Snow Band File
			area_frac_band[cell][band_idx] = band_areas[cell][band_idx] / total_cell_area
			# This will be used to update the Vegetation Parameter File
			area_frac_glacier[cell][band_idx] = glacier_areas[cell][band_idx] / total_cell_area
	return area_frac_band, area_frac_glacier


if __name__ == '__main__':

	print '\n\nVIC + RGM ... together at last!'
	# Get all initial VIC global parameters
	global_parms = get_global_parms(vic_global_file)
	# Get entire time range of coupled VIC-RGM run from the initial VIC global file
	start_year = int(global_parms['STARTYEAR'][0][0])
	end_year = int(global_parms['ENDYEAR'][0][0])
	# Initial VIC output state filename prefix is determined by STATENAME in the global file
	state_filename_prefix = str(global_parms['STATENAME'][0][0])

	# Get VIC vegetation parameters and grid cell IDs from initial Vegetation Parameter File
	veg_parm_file = global_parms['VEGPARAM'][0][0]
	cell_ids, num_veg_tiles, veg_parms = get_veg_parms(veg_parm_file)
	# Get VIC snow/elevation band parameters from initial Snow Band File
	num_snow_bands = int(global_parms['SNOW_BAND'][0][0])
	snb_file = global_parms['SNOW_BAND'][0][1]
	snb_parms = get_snb_parms(snb_file, num_snow_bands)
	# Get list of elevation bands for each VIC grid cell
	band_map = get_bands(snb_parms, band_size)

	# The RGM will always output a DEM file of the same name (if running RGM for a single year at a time)
	rgm_surf_dem_out_file = temp_files_path + 's_out_00001.grd'

	# Open and read VIC-grid-to-RGM-pixel mapping file
	# pixel_to_cell_map is a list of dimensions num_rows_rpg x num_cols_rpg, each element containing a VIC grid cell ID
	pixel_to_cell_map, num_rows_rpg, num_cols_rpg = get_rgm_pixel_mapping(pixel_cell_map_file)

	# Open and read the provided Bed Digital Elevation Map (BDEM) file into a 2D bdem_grid
	# Also get RGM pixel area in meters, and min/max x and y coordinates 
	bed_dem, rgm_pixel_area, rgm_xmin, rgm_xmax, rgm_ymin, rgm_ymax = read_asc_dem_file(bed_dem_file)


	# NOTE: it appears the .gsa DEM from Markus does not conform to ARC/ASCII grid format (rows are clipped at 10 elements), 
	# which makes RGM puke.  Once this is corrected, we can eliminate the need for the calls to write_grid_to_file() to 
	# write out the Bed and Surface DEMs to a format readable by RGM.  Ideally, we probably only want to have to read .gsa files, 
	# and remove support / need for .asc files entirely
	rgm_bed_dem_file = input_files_path + 'peyto_bed_dem.gsa'  # this file will now point at the ARC/ASCII version, produced by read_pcic_dem_file, of the Bed DEM Markus provided (temporary workaround)
	write_grid_to_gsa_file(bed_dem, rgm_bed_dem_file)
	rgm_surf_dem_file = input_files_path + 'peyto_surf_dem.gsa'
	write_grid_to_gsa_file(bed_dem, rgm_surf_dem_file) # using the bed DEM as the initial surface DEM for testing


	for year in range(start_year, end_year):
		# Call initial VIC run starting at the first year in the VIC global parameters file
		print '\nRunning year: {}'.format(year)

		# 1. Write / Update temporary Global Parameters File, temp_gpf
		# If year > start_year, introduce INIT_STATE line with most current state_file (does not exist in the first read-in of global_parms).
		# Overwrite VEGPARAM parameter with temp_vpf, and SNOW_BAND with temp_snb
		state_file = state_filename_prefix + "_" + str(year) + str(global_parms['STATEMONTH'][0][0]) + str(global_parms['STATEDAY'][0][0])
		update_global_parms(year > start_year)
		temp_gpf = write_global_parms_file()
		print 'invoking VIC with global parameter file {}'.format(temp_gpf)

		# 2. Run VIC for a year.  This will save VIC model state at the end of the year, along with a Glacier Mass Balance (GMB) polynomial for each cell
		subprocess.check_call([vic_full_path, "-g", temp_gpf], shell=False, stderr=subprocess.STDOUT)

		# 3. Open VIC NetCDF state file and get the most recent GMB polynomial for each grid cell being modeled
		print 'opening state file {}'.format(state_file)
		state = h5py.File(state_file, 'r+')
		gmb_polys = get_mass_balance_polynomials(state, state_file)
			
		# 4. Translate mass balances using grid cell GMB polynomials and current veg_parm_file into a 2D RGM mass balance grid (MBG)
		mass_balance_grid = mass_balances_to_rgm_grid(gmb_polys, pixel_to_cell_map)
		# write Mass Balance Grid to ASCII file to direct the RGM to use as input
		mbg_file = temp_files_path + 'mass_balance_grid_' + str(year) + '.gsa'
		write_grid_to_gsa_file(mass_balance_grid, mbg_file)

		# 5. Run RGM for one year, passing MBG, BDEM, SDEM
		#subprocess.check_call([rgm_full_path, "-p", rgm_params_file, "-b", bed_dem_file, "-d", sdem_file, "-m", mbg_file, "-o", temp_files_path, "-s", "0", "-e", "0" ], shell=False, stderr=subprocess.STDOUT)
		subprocess.check_call([rgm_full_path, "-p", rgm_params_file, "-b", rgm_bed_dem_file, "-d", rgm_surf_dem_in_file, "-m", mbg_file, "-o", temp_files_path, "-s", "0", "-e", "0" ], shell=False, stderr=subprocess.STDOUT)
		# remove temporary files if not saving for offline inspection
		if not output_trace_files:
			os.remove(mbg_file)
			os.remove(rgm_surf_dem_file)

		# 6. Read in new Surface DEM file from RGM output
		rgm_surf_dem_out, sdem_pixel_area, sdem_xmin, sdem_xmax, sdem_ymin, sdem_ymax = read_gsa_dem_file(rgm_surf_dem_out_file)
		temp_surf_dem_file = temp_files_path + 'rgm_surf_dem_out_' + str(year) + '.gsa'
		os.rename(rgm_surf_dem_out_file, temp_surf_dem_file)
		# this will be fed back into RGM on next time step
		rgm_surf_dem_in_file = temp_surf_dem_file

		# 7. Update glacier mask
		glacier_mask = update_glacier_mask(rgm_surf_dem_out, bed_dem)
		if output_trace_files:
			glacier_mask_file = temp_files_path + 'glacier_mask_' + str(year) + '.gsa'
			write_grid_to_gsa_file(glacier_mask, glacier_mask_file)
		
		# 8. Update areas of each elevation band in each VIC grid cell, and calculate area fractions
		area_frac_bands, area_frac_glacier = update_band_areas()

		# 9. Update vegetation parameters and write to new temporary file temp_vpf
		update_veg_parms()
		temp_vpf = write_veg_parms_file()

		# 10. Update snow band parameters and write to new temporary file temp_snb
		update_snb_parms()
		temp_snb = write_snb_parms_file()

		# 11 Update HRUs in VIC state file 
			# don't forget to close the state file
		