#!/usr/bin/env python

""" This script orchestrates a coupled Variable Infiltration Capacity (VIC) and Regional Glacier Model (RGM) run. """

import sys
import argparse
import subprocess
import numpy as np
import h5py


vic_full_path = '/home/mfischer/code/vic/vicNl'
#initial_vic_global_file = '/home/mfischer/vic_dev/input/place/glb_base_PLACE_19601995_VIC4.1.2_outNETCDF_initial.txt' 
#vic_global_file = '/home/mfischer/vic_dev/input/place/glb_base_PLACE_19601995_VIC4.1.2_outNETCDF.txt' 
#vpf_full_path = '/home/mfischer/vic_dev/input/place/vpf_place_100m.txt' 
# path to Glacier Mass Balance files
gmb_path = '/home/mfischer/vic_dev/out/testing/gmb_files/'


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s]n' % message)
        self.print_help()
        sys.exit(2)

def get_global_parameters(global_parm_file):
	""" Parses the VIC global parameters file created by the user that will be needed to conduct the VIC-RGM run """
	# Important parms: STARTYEAR, ENDYEAR, GLACIER_ID, GLACIER_ACCUM_START_YEAR, GLACIER_ACCUM_INTERVAL, STATENAME, VEGPARAM, SOIL, SNOW_BAND, RESULT_DIR
	global_parms = {}
	with open(global_parm_file, 'r') as f:
		for line in iter(f):
			#print 'line: x{}x'.format(line)
			if not line.isspace() and line[0] is not '#':
				columns = line.split()
				#print 'columns: {}'.format(columns)
				num_columns = len(columns)
				parm_name = columns[0]
				try:
					if global_parms[parm_name]: # if we've already read one or more entries of this parm_name
#						print 'parm {} exists already'.format(parm_name)
						for i in range(1, num_columns):
							temp_list.append(columns[i])
						global_parms[parm_name].append(temp_list)
				except:
					global_parms[parm_name] = []
					temp_list = []
					for i in range(1, num_columns):
						temp_list.append(columns[i])
					global_parms[parm_name].append(temp_list)
#				print 'global_parms[{}]: {}'.format(parm_name,global_parms[parm_name])
		f.close()
	return global_parms

def get_grid_cell_ids(veg_parm_file):
	""" Parses the vegetation parameter file for the grid cell IDs that will be used in the simulation """
	cell_ids = []
	with open(veg_parm_file, 'r') as f:
		for line in iter(f):
			columns = line.split()
			if len(columns) == 2:
				cell_ids.append(columns[0])
		f.close()
	return cell_ids

def get_rgm_to_grid_cell_mapping(rgm_vic_map_file):
	""" Parses the RGM-to-VIC pixel-to-grid cell map file and initialises a 2D map of valid RGM pixels to use as one input to RGM """
	pass

def get_gmb_polynomials(gmb_file, year, cell_ids):
	""" Opens glacier mass balance file and retrieves the most recent GMB polynomial for each grid cell in cell_ids """
		# GMB file format to be written out by VIC:
		# <year> <month> <day> <cell ID> <b0> <b1> <b2>
		pass

def mass_balances_to_rgm_pixels(gmb_polys, rgm_pixel_map):
	""" Translate mass balances from grid cell GMB polynomials to 2D RGM pixel grid to use as one input to RGM """
		pass

def write_mbg_to_file(mass_balance_grid):
	""" Writes the 2D Mass Balance Grid to ASCII file in the input format expected by the RGM, and returns a filename """
		pass

if __name__ == '__main__':
	
	parser = MyParser()
	parser.add_argument('-g', action="store", dest="vic_global_file", type=str, help = 'file name and path of the VIC global parameters file')
	parser.add_argument('-vpf', action="store", dest="veg_parm_file", type=str, help = 'file name and path of the initial Vegetation Parameter File (VPF)')
	parser.add_argument('-sdem', action="store", dest="surf_dem_file", type=str, help = 'file name and path of the initial Surface Digital Elevation Model (SDEM) file')
	parser.add_argument('-bdem', action="store", dest="bed_dem_file", type=str, help = 'file name and path of the Bed Digital Elevation Model (BDEM) file')

	if len(sys.argv) == 1:
    	parser.print_help()
    	sys.exit(1)
	options = parser.parse_args()
	vic_global_file = options.vic_global_file
	initial_veg_parm_file = options.veg_parm_file
	initial_surf_dem_file = options.surf_dem_file
	bed_dem_file = options.bed_dem_file

	# Get all VIC global parameters
	global_parms = get_global_parameters(vic_global_file)

	# Get all grid cell ID numbers for this simulation from the initial Vegetation Parameter File (iVPF)
	cell_ids = get_grid_cell_ids(global_parms['VEGPARAM'])

	# Open and read RGM-VIC pixel-cell mapping file
	rgm_pixel_map = get_rgm_to_grid_cell_mapping()

	# Open and read Bed Digital Elevation Map (BDEM) file into 2D bdem_grid
	bdem_grid = get_bdem(bed_dem_file)

	start_year = global_parms['STARTYEAR']
	end_year = global_parms['ENDYEAR']
	global_file = vic_global_file
	veg_file = initial_veg_parm_file
	sdem_file = initial_surf_dem_file

	for year in range(start_year, end_year):
		# Call initial VIC run starting at year 0. 
		if year > start_year:
			global_file = temp_gpf
			veg_file = temp_vpf
			sdem_file = temp_sdem

		# Run VIC for a year.  This will do the following:
		# 1) create, if needed, and write a polynomial line to the Glacier Mass Balance (GMB) file
		# 2) save VIC state at the end of the year
		subprocess.check_call([vic_full_path, "-g", global_file], shell=False, stderr=subprocess.STDOUT)

		# 3. Open VIC state file and get the most recent GMB polynomial for each grid cell being modeled
		state_filename = global_parms['STATENAME'] + "_" + year + global_parms['STATEMONTH'] + global_parms['STATEDAY']
		state = h5py.open(state_filename)
		gmb_polys = state['GLAC_MASS_BALANCE_EQUATION']
			
		# 4. Translate mass balances using grid cell GMB polynomials and current veg_parm_file into a 2D RGM mass balance grid (MBG)
		mass_balance_grid = mass_balances_to_rgm_pixels(gmb_polys, veg_parm_file, rgm_pixel_map)
		# write MBG to ASCII file to direct the RGM to use as input
		mass_balance_file = write_mbg_to_file(mass_balance_grid)

		# 5. Run RGM for one year, passing MBG, BDEM, SDEM
		subprocess.check_call([rgm_full_path, "-p", rgm_params_file, "-b", bed_dem_file, "-d", sdem_file, "-m" mass_balance_file, "-s", 0, "-e", 0 ], shell=False, stderr=subprocess.STDOUT)
	
		# 6. Read in new Surface DEM file from RGM output
		sdem_file = get_temp_sdem_filename()
		temp_sdem = read_sdem_file()

		# 7. Update glacier mask
		glacier_mask_grid = update_glacier_mask(temp_sdem, bdem_grid)
		
		# 8.1 Update HRUs in VIC state file 

		# 8.2 Write / Update temp_vpf

		# 9. Write / Update temp_gpf