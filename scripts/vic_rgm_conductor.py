#!/usr/bin/env python

""" This script orchestrates a coupled Variable Infiltration Capacity (VIC) 
  and Regional Glacier Model (RGM) run. 
"""

import argparse
import csv
import collections
from collections import OrderedDict
import os
import shutil
import subprocess
import sys
from warnings import warn

import numpy as np
import h5py
from dateutil.relativedelta import relativedelta

from conductor.io import get_rgm_pixel_mapping, read_gsa_headers,\
  write_grid_to_gsa_file, update_glacier_mask, read_state
from conductor.cells import Band, HydroResponseUnit
from conductor.snbparams import load_snb_parms, save_snb_parms
from conductor.vegparams import load_veg_parms, save_veg_parms
from conductor.vic_globals import Global

# These should become command line parameters at some point
vic_full_path = '/home/mfischer/code/vic/vicNl'
rgm_full_path = '/home/mfischer/code/rgm/rgm'
temp_files_path = '/home/mfischer/vic_dev/out/testing/hydrocon_temp/'
# set it as default = os.env(tmp)

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
  parser.add_argument('--g', action="store", dest="vic_global_file", type=str,\
      help = 'file name and path of the VIC global parameters file')
  parser.add_argument('--rgm-params', action="store", dest="rgm_params_file",\
    type=str, help = 'file name and path of the Regional Glacier Model (RGM)\
    parameters file')
  parser.add_argument('--sdem', action="store", dest="surf_dem_file", type=str,\
    help = 'file name and path of the initial Surface Digital Elevation Model \
    (SDEM) file (GSA format)')
  parser.add_argument('--bdem', action="store", dest="bed_dem_file", type=str,\
    help = 'file name and path of the Bed Digital Elevation Model (BDEM) file \
    (GSA format)')
  parser.add_argument('--pixel-map', action="store", dest="pixel_cell_map_file",\
    type=str, help = 'file name and path of the RGM Pixel to VIC Grid Cell \
    mapping file')
  parser.add_argument('--glacier-mask', action="store",\
    dest="init_glacier_mask_file", type=str, help = 'file name and path of the \
    file containing the initial glacier mask (GSA format)')
  parser.add_argument('--trace-files', action="store_true", default=False,\
    dest="trace_files", help = 'write out persistent GSA format surface DEMs \
    and glacier masks, and 2D mass balance grid files, on each time step for \
    offline inspection')
  parser.add_argument('--open-ground-root-zone', action="store",\
    dest="open_ground_root_zone_file", type=str, default=None,\
    help = 'file name and path of one-line text file containing 6 custom root \
    parameters for the bare soil vegetation type / HRU (same format as a \
      vegetation tile line in the vegetation parameters file). \
      Default: 0.10  1.00  0.10  0.00  0.10  0.00')
  parser.add_argument('--glacier-root-zone', action="store",\
    dest="glacier_root_zone_file", type=str, default=None,\
    help = 'file name and path of one-line text file containing 6 custom root \
    parameters for the glacier vegetation type / HRU (same format as a \
      vegetation tile line in the vegetation parameters file). \
      Default: 0.10  1.00  0.10  0.00  0.10  0.00')
  parser.add_argument('--band-size', action="store", dest="band_size", \
    type=int, default=100,\
    help = 'vertical size of VIC elevation bands in metres (default = 100m)')

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

  if open_ground_root_zone_file:
    with open(open_ground_root_zone_file, 'r') as f:
      line = f.readline()
      open_ground_root_zone_parms = [float(x) for x in line.split()]
      if len(open_ground_root_zone_parms) != 6:
        print('Open ground root zone parameters file is malformed. Expected \
          6 space-separated numeric values on a single line. Exiting.\n')
        sys.exit(0)
  else:
    open_ground_root_zone_parms = None

  if glacier_root_zone_file:
    with open(glacier_root_zone_file, 'r') as f:
      line = f.readline()
      glacier_root_zone_parms = [float(x) for x in line.split()]
      if len(glacier_root_zone_parms) != 6:
        print('Glacier root zone parameters file is malformed. Expected \
          6 space-separated numeric values on a single line. Exiting.\n')
        sys.exit(0)
  else:
    glacier_root_zone_parms = None

  return vic_global_file, rgm_params_file, surf_dem_in_file, bed_dem_file, \
    pixel_cell_map_file, init_glacier_mask_file, output_trace_files, \
    glacier_root_zone_parms, open_ground_root_zone_parms, band_size

def mass_balances_to_rgm_grid(gmb_polys, pixel_to_cell_map, num_rows_dem,\
  num_cols_dem, cell_ids):
  """ Translate mass balances from grid cell GMB polynomials to 2D RGM pixel \
    grid to use as one of the inputs to RGM
  """
  mass_balance_grid = [[0 for x in range(num_cols_dem)]\
    for x in range(num_rows_dem)]
  try:
    for row in range(num_rows_dem):
      for col in range(num_cols_dem):
        pixel = pixel_to_cell_map[row][col]
        # band = pixel[0]
        cell_id = pixel[0]
        # read most recent median elevation of this pixel
        median_elev = pixel[1]
        # only grab pixels that fall within a VIC cell
        if cell_id != 'NA':
          # check that the cell_id agrees with reading from the veg_parm_file
          if cell_id not in cell_ids:
            print('mass_balances_to_rgm_grid: Cell ID {} was not found in the \
              list of VIC cell IDs read from the vegetation parameters file. \
              Exiting.\n'.format(cell_id))
            sys.exit(0)
          mass_balance_grid[row][col] = gmb_polys[cell_id][0] + median_elev\
            * (gmb_polys[cell_id][1] + median_elev * gmb_polys[cell_id][2])
  except:
    print('mass_balances_to_rgm_grid: Error while processing pixel {} \
      (row {} column {})'.format(pixel, row, col))
  return mass_balance_grid

def run_ranges(startdate, enddate, glacier_start):
  """Generator which yields date ranges (a 2-tuple) that represent times at
    which to begin and end a VIC run.
    startdate and enddate are the overall boundaries of the simulation.
    glacier_start is a date some time after startdate that *should* be aligned
    with the water year (this code mostly assumes that that is the case).
    Arguments (1950/01/01, 1995/12/31, 1955/10/01) would yield a sequence
    as follows:
    (1950/01/01, 1956/09/30)
    (1956/10/01, 1957/09/30)
    (1957/10/01, 1958/09/30)
    (1958/10/01, 1959/09/30)
    ...
    (1994/10/01, 1995/09/30)
    Note that the sequence will and should end before the specified end date,
    aligned with the water year.
  """
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
    glacier_root_zone_parms, open_ground_root_zone_parms, band_size\
    = parse_input_parms()

  # Get all initial VIC global parameters from the global parameter file
  with open(vic_global_file, 'r') as f:
    global_parms = Global(f)

  assert global_parms.state_format == 'NETCDF', \
    "{} only supports NetCDF input statefile input as opposed "\
    "to the specified {}. Please change change this in your "\
    "global file {}".format(__name__, global_parms.state_format,\
    vic_global_file)

  # Initial VIC output state filename prefix is determined by STATENAME
  # in the global file
  state_filename_prefix = global_parms.statename

  # Apply custom glacier_id and open_ground_id Band attributes, if provided
  if global_parms.glacier_id is None:
    print('No value for GLACIER_ID was provided in the VIC global file. \
      Assuming default value of {}.'.format(Band.glacier_id))
  else:
    Band.glacier_id = global_parms.glacier_id
  # FIXME: reinstate the following commented-out code once OPEN_GROUND_ID
  # is supported in VIC
  # Numeric code indicating an open ground vegetation tile (HRU)
  # if global_parms.open_ground_id is None:
  #     print('No value for OPEN_GROUND_ID was provided in the VIC global file.
  # Assuming default value of {}.'.format(Band.open_ground_id))
  # else:
  #     Band.open_ground_id = global_parms.open_ground_id
  if band_size:
    Band.band_size = band_size

  # Load parameters from Snow Band Parameters File
  num_snow_bands, snb_file = global_parms.snow_band.split()
  num_snow_bands = int(num_snow_bands)
  elevation_cell_dict = load_snb_parms(snb_file, num_snow_bands)

  # Load vegetation parameters from initial Vegetation Parameter File
  hru_cell_dict = vegparams.load_veg_parms(global_parms.vegparam)

  # Apply custom HRU root_zone_parms attributes, if provided
  if glacier_root_zone_parms or open_ground_root_zone_parms:
    cells.apply_custom_root_zone_parms(hru_cell_dict, glacier_root_zone_parms,\
      open_ground_root_zone_parms)
    Band.glacier_root_zone_parms = glacier_root_zone_parms
    Band.open_ground_root_zone_parms = open_ground_root_zone_parms

  # Create Ordered dictionary of Cell objects by merging info gathered from
  # Snow Band and Vegetation Parameter files and custom parameters
  cells = merge_cell_input(hru_cell_dict, elevation_cell_dict)

  # TODO: Do a sanity check to make sure band area fractions in Snow Band
  # Parameters file add up to sum of HRU area fractions in Vegetation
  # Parameter File for each cell?
  #assert (area_fracs == [sums of HRU area fracs for all bands])
  
  # The RGM will always output a DEM file of the same name (if running RGM for
  # a single year at a time)
  rgm_surf_dem_out_file = temp_files_path + 's_out_00001.grd'

  # Open and read VIC-grid-to-RGM-pixel mapping file
  cellid_map, elevation_map, cell_areas, num_cols_dem, num_rows_dem\
    = get_rgm_pixel_mapping(pixel_cell_map_file)

  # Get DEM xmin, xmax, ymin, ymax metadata of Bed DEM and check file header
  # validity     
  dem_xmin, dem_xmax, dem_ymin, dem_ymax, num_rows, num_cols\
    = read_gsa_headers(bed_dem_file)
  # Verify that number of columns & rows agree with what's stated in the
  # pixel_to_cell_map_file
  assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch \
    of stated dimension(s) between Bed DEM in {} (num rows: {}, num columns: \
    {}) and the VIC-grid-to-RGM-pixel map in {} (num rows: {}, num columns: \
    {}). Exiting.\n'.format(bed_dem_file, num_rows, num_cols,\
    pixel_cell_map_file, num_rows_dem, num_cols_dem)

  # Read in the provided Bed Digital Elevation Map (BDEM) file to 2D bed_dem array
  bed_dem = np.loadtxt(bed_dem_file, skiprows=5)

  # Check header validity of Surface DEM file
  _, _, _, _, num_rows, num_cols = read_gsa_headers(surf_dem_in_file)
  # Verify number of columns & rows agree with what's stated in the
  # pixel_to_cell_map_file
  assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch \
    of stated dimension(s) between Surface DEM in {} (num rows: {}, \
    num columns: {}) and the VIC-grid-to-RGM-pixel map in {} (num rows: {}, \
    num columns: {}). Exiting.\n'.format(surf_dem_in_file, num_rows, num_cols,\
    pixel_cell_map_file, num_rows_dem, num_cols_dem)

  # Read in the provided Surface Digital Elevation Map (SDEM) file to 2D surf_dem array
  surf_dem_initial = np.loadtxt(surf_dem_in_file, skiprows=5)
  # Check agreement between elevation map from VIC-grid-to-RGM-pixel file and
  # the initial Surface DEM
  assert (np.equal(elevation_map, surf_dem_initial)), 'Values mismatch \
    between provided initial Surface DEM file (num rows:{}, num columns: \
    {}) and VIC-grid-to-RGM-pixel file (num rows: {}, num columns: {}). \
    Exiting.\n'.format(num_rows, num_cols, num_rows_dem, num_cols_dem)

  # Check header validity of initial Glacier Mask file
  _, _, _, _, num_rows, num_cols = read_gsa_headers(init_glacier_mask_file)
  # Verify number of columns & rows agree with what's stated in the 
  # pixel_to_cell_map_file
  assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch \
  of stated dimension(s) between Glacier Mask in {} (num rows: {}, \
  num columns: {}) and the VIC-grid-to-RGM-pixel map in {} (num rows: {}, num \
  columns: {}). Exiting.\n'.format(init_glacier_mask_file, num_rows, num_cols,\
  pixel_cell_map_file, num_rows_dem, num_cols_dem)
  # Read in the provided initial glacier mask file to 2D glacier_mask array 
  glacier_mask = np.loadtxt(init_glacier_mask_file, skiprows=5)

  # Apply the initial glacier mask and modify the band and glacier area
  # fractions accordingly
  update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,
    surf_dem_initial, num_rows_dem, num_cols_dem, glacier_mask)

  # Write temporary snow band and vegetation parameter files for first VIC run
  temp_snb = temp_files_path + 'snb_temp_'\
    + global_parms.startdate.isoformat() + '.txt'
  snbparams.save_snb_parms(cells, temp_snb, band_map)
  temp_vpf = temp_files_path + 'vpf_temp_'\
    + global_parms.startdate.isoformat() + '.txt'
  vegparams.save_veg_parms(cells, temp_vpf)

# (initialisation done)

#### Run the coupled VIC-RGM model for the time range specified in the VIC
  # global parameters file
  time_iterator = run_ranges(global_parms.startdate,
                 global_parms.enddate,
                 global_parms.glacier_accum_startdate)
  for start, end in time_iterator:
    # 1. Write / Update temporary Global Parameters File, temp_gpf
    temp_gpf = temp_files_path + 'gpf_temp_{}.txt'.format(start.isoformat())
    print('\nRunning VIC from {} to {}, using global parameter file {}'.format(start, end, temp_gpf))
    # set global parameters for this VIC run
    global_parms.vegparm = temp_vpf
    global_parms.snow_band = '{} {}'.format(num_snow_bands, temp_snb)
    global_parms.startdate = start
    global_parms.enddate = end
    global_parms.statedate = end

    global_parms.write(temp_gpf)

    # 2. Run VIC for a year.  This will save VIC model state at the end of the
    # year, along with a Glacier Mass Balance (GMB) polynomial for each cell
    subprocess.check_call([vic_full_path, "-g", temp_gpf], shell=False,\
      stderr=subprocess.STDOUT)

    # 3. Open VIC NetCDF state file and load the most recent set of state
    # variable values for all grid cells being modeled
    state_file = state_filename_prefix + "_" + start.isoformat() ## FIXME
    print('opening VIC state file {}'.format(state_file))
    # leave the state file open for modification later
    state = h5py.File(state_file, 'r+')
    # read new states of all cells
    read_state(state, cells)
    if output_trace_files: # keep a copy of the original state file from VIC
        shutil.copy(state_file, state_file + '_orig')
    # pull Glacier Mass Balance polynomials out of cell states
    gmb_polys = {}
    for cell_id in cells:
        # leave off the "fit error" term at the end of the GMB polynomial
        # TODO: generalize this to handle variable length GMB polynomials
        gmb_polys[cell_id] = cells[cell_id].cell_state.variables['GLAC_MASS_BALANCE_INFO'][0:3]

    # 4. Translate mass balances using grid cell GMB polynomials and current
    # veg_parm_file into a 2D RGM mass balance grid (MBG)
    mass_balance_grid = mass_balances_to_rgm_grid(gmb_polys, pixel_to_cell_map,\
      num_rows_dem, num_cols_dem, cell_ids)
    # write Mass Balance Grid to ASCII file to direct the RGM to use as input
    mbg_file = temp_files_path + 'mass_balance_grid_' + start.isoformat()\
      + '.gsa'
    write_grid_to_gsa_file(mass_balance_grid, mbg_file, num_cols_dem,\
      num_rows_dem, dem_xmin, dem_xmax, dem_ymin, dem_ymax)

    # 5. Run RGM for one year, passing it the MBG, BDEM, SDEM
    subprocess.check_call([rgm_full_path, "-p", rgm_params_file, "-b",\
      bed_dem_file, "-d", surf_dem_in_file, "-m", mbg_file, "-o",\
      temp_files_path, "-s", "0", "-e", "0" ], shell=False,\
      stderr=subprocess.STDOUT)
    # remove temporary files if not saving for offline inspection
    if not output_trace_files:
      os.remove(mbg_file)
      os.remove(rgm_surf_dem_file)

    # 6. Read in new Surface DEM file from RGM output
    rgm_surf_dem_out = np.loadtxt(rgm_surf_dem_out_file, skiprows=5)
    temp_surf_dem_file = temp_files_path + 'rgm_surf_dem_out_'\
      + start.isoformat() + '.gsa'
    os.rename(rgm_surf_dem_out_file, temp_surf_dem_file)
    # this will be fed back into RGM on next time step:
    surf_dem_in_file = temp_surf_dem_file

    # 7. Update glacier mask
    glacier_mask = update_glacier_mask(rgm_surf_dem_out, bed_dem,\
      num_rows_dem, num_cols_dem)
    if output_trace_files:
      glacier_mask_file = temp_files_path + 'glacier_mask_'\
        + start.isoformat() + '.gsa'
      write_grid_to_gsa_file(glacier_mask, glacier_mask_file)
    
    # 8. Update areas of each elevation band in each VIC grid cell, and update
    # snow band and vegetation parameters
    update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands, \
      rgm_surf_dem_out, num_rows_dem, num_cols_dem, glacier_mask)
    temp_snb = temp_files_path + 'snb_temp_' + start.isoformat() + '.txt'
    snbparams.save_snb_parms(cells, temp_snb, band_map)
    temp_vpf = temp_files_path + 'vpf_temp_' + start.isoformat() + '.txt'
    vegparams.save_veg_parms(cells, temp_vpf)

    # 11 Calculate changes to HRU state variables and modify the VIC state file
    #update_state(cells)
    #write_state(cells, state)
    state.close()

    # Get ready for the next loop
    global_parms.init_state = "{}_{}.txt".format(global_parms.statename,\
      end.strftime("%Y%m%d"))
    global_parms.statedate = end

# Main program invocation.
if __name__ == '__main__':
  main()
