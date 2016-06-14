#!/usr/bin/env python

""" This script orchestrates a coupled Variable Infiltration Capacity (VIC) 
  and Regional Glacier Model (RGM) run. 
"""

import argparse
import collections
from collections import OrderedDict
import os
import shutil
import subprocess
import sys
from warnings import warn
import logging

import numpy as np
import netCDF4
from dateutil.relativedelta import relativedelta

from conductor.io import get_rgm_pixel_mapping, read_gsa_headers,\
  write_grid_to_gsa_file, update_glacier_mask, read_state, write_state
from conductor.cells import Cell, Band, HydroResponseUnit, merge_cell_input,\
  update_area_fracs
from conductor.snbparams import load_snb_parms, save_snb_parms
from conductor.vegparams import load_veg_parms, save_veg_parms
from conductor.vic_globals import Global

# These should become command line parameters at some point
vic_full_path = '/home/mfischer/code/vic/vicNl'
rgm_full_path = '/home/mfischer/code/rgm/rgm'
output_path = '/home/mfischer/vic_dev/out/testing/'
temp_files_path = output_path + 'hydrocon_temp/'
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
  parser.add_argument('--g', action='store', dest='vic_global_file', type=str,\
      help='file name and path of the VIC global parameters file')
  parser.add_argument('--rgm-params', action='store', dest='rgm_params_file',\
    type=str, help='file name and path of the Regional Glacier Model (RGM)\
    parameters file')
  parser.add_argument('--sdem', action='store', dest='surf_dem_file', type=str,\
    help='file name and path of the initial Surface Digital Elevation Model \
    (SDEM) file (GSA format)')
  parser.add_argument('--bdem', action='store', dest='bed_dem_file', type=str,\
    help='file name and path of the Bed Digital Elevation Model (BDEM) file \
    (GSA format)')
  parser.add_argument('--pixel-map', action='store', dest='pixel_cell_map_file',\
    type=str, help='file name and path of the RGM Pixel to VIC Grid Cell \
    mapping file')
  parser.add_argument('--glacier-mask', action='store',\
    dest='init_glacier_mask_file', type=str, help='file name and path of the \
    file containing the initial glacier mask (GSA format)')
  parser.add_argument('--trace-files', action='store_true', default=False,\
    dest='trace_files', help='write out persistent GSA format surface DEMs \
    and glacier masks, and 2D mass balance grid files, on each time step for \
    offline inspection')
  parser.add_argument('--open-ground-root-zone', action='store',\
    dest='open_ground_root_zone_file', type=str, default=None,\
    help='file name and path of one-line text file containing 6 custom root \
    parameters for the bare soil vegetation type / HRU (same format as a \
      vegetation tile line in the vegetation parameters file). \
      Default: 0.10  1.00  0.10  0.00  0.10  0.00')
  parser.add_argument('--glacier-root-zone', action='store',\
    dest='glacier_root_zone_file', type=str, default=None,\
    help='file name and path of one-line text file containing 6 custom root \
    parameters for the glacier vegetation type / HRU (same format as a \
      vegetation tile line in the vegetation parameters file). \
      Default: 0.10  1.00  0.10  0.00  0.10  0.00')
  parser.add_argument('--band-size', action='store', dest='band_size',\
    type=int, default=100,\
    help='vertical size of VIC elevation bands in metres (default = 100m)')
  parser.add_argument('--loglevel', action='store', dest='loglevel', type=str,\
    default='INFO', help='the logging verbosity level to be written to the \
      log file. Options are: DEBUG, INFO, WARNING, ERROR.')

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
  loglevel = options.loglevel

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
    glacier_root_zone_parms, open_ground_root_zone_parms, band_size, loglevel

def mass_balances_to_rgm_grid(gmb_polys, cell_id_map, dem, num_rows_dem,\
  num_cols_dem, valid_cell_ids):
  """ Translate mass balances from grid cell GMB polynomials to 2D RGM pixel \
    grid to use as one of the inputs to RGM
  """
  mass_balance_grid = [[0 for x in range(num_cols_dem)]\
    for x in range(num_rows_dem)]
  try:
    for row in range(num_rows_dem):
      for col in range(num_cols_dem):
        if not np.isnan(cell_id_map[row][col]):
        # only grab elevation for pixels that fall within a VIC cell
          cell_id = str(int(cell_id_map[row][col]))
          # check that the cell_id agrees with reading from the veg_parm_file
          if cell_id not in valid_cell_ids:
            print('mass_balances_to_rgm_grid: Cell ID {} was not found in the \
              list of VIC cell IDs read from the vegetation parameters file. \
              Exiting.\n'.format(cell_id))
            sys.exit(0)
          # read most recent median elevation of this pixel
          median_elev = dem[row][col]
          mass_balance_grid[row][col] = gmb_polys[cell_id][0] + median_elev \
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
  vic_global_file, rgm_params_file, surf_dem_in_file, bed_dem_file,\
    pixel_cell_map_file, init_glacier_mask_file, output_trace_files,\
    glacier_root_zone_parms, open_ground_root_zone_parms, band_size, loglevel\
    = parse_input_parms()

  # Set up logging
  numeric_loglevel = getattr(logging, loglevel.upper())
  logging.basicConfig(filename=output_path+'hydrocon.log',\
    level=numeric_loglevel, format='%(levelname)s %(asctime)s %(message)s')
  logging.info('------- VIC-RGM Hydro-Conductor Startup -------')

  # Get all initial VIC global parameters from the global parameter file
  logging.info('Loading initial VIC global parameters from %s', vic_global_file)
  with open(vic_global_file, 'r') as f:
    global_parms = Global(f)

  assert global_parms.state_format == 'NETCDF',\
    "VIC only supports NETCDF input statefile input format."\
    "(STATE_FORMAT in global file is currently set as {})"\
    .format(global_parms.state_format)

  # Apply custom glacier_id and open_ground_id Band attributes, if provided
  if global_parms.glacier_id is None:
    logging.info('No value for GLACIER_ID was provided in the VIC global file. \
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
  logging.info('Loading initial VIC snow band parameters from %s', snb_file)
  elevation_cell_dict = load_snb_parms(snb_file, num_snow_bands)

  # Load vegetation parameters from initial Vegetation Parameter File
  logging.info('Loading initial VIC vegetation parameters from %s',\
    global_parms.vegparam)
  hru_cell_dict = load_veg_parms(global_parms.vegparam)

  # Apply custom HRU root_zone_parms attributes, if provided
  if glacier_root_zone_parms or open_ground_root_zone_parms:
    cells.apply_custom_root_zone_parms(hru_cell_dict, glacier_root_zone_parms,\
      open_ground_root_zone_parms)
    Band.glacier_root_zone_parms = glacier_root_zone_parms
    Band.open_ground_root_zone_parms = open_ground_root_zone_parms

  # TODO: Do a sanity check to make sure band area fractions in Snow Band
  # Parameters file add up to sum of HRU area fractions in Vegetation
  # Parameter File for each cell?
  #assert (area_fracs == [sums of HRU area fracs for all bands])

  # Create Ordered dictionary of Cell objects by merging info gathered from
  # Snow Band and Vegetation Parameter files and custom parameters
  cells = merge_cell_input(hru_cell_dict, elevation_cell_dict)

  # Open and read VIC-grid-to-RGM-pixel mapping file
  logging.info('Loading VIC-grid-to-RGM-pixel mapping from %s',\
    pixel_cell_map_file)
  cell_id_map, elevation_map, cell_areas, num_cols_dem, num_rows_dem\
    = get_rgm_pixel_mapping(pixel_cell_map_file)

  # Get DEM xmin, xmax, ymin, ymax metadata of Bed DEM and check file header
  # validity     
  dem_xmin, dem_xmax, dem_ymin, dem_ymax, num_rows, num_cols\
    = read_gsa_headers(bed_dem_file)
  # Verify that number of columns & rows agree with what's stated in the
  # pixel_cell_map_file
  assert (num_cols == num_cols_dem) and (num_rows == num_rows_dem),'Mismatch \
  of stated dimension(s) between Bed DEM in {} (num rows: {}, num columns: \
  {}) and the VIC-grid-to-RGM-pixel map in {} (num rows: {}, num columns: \
  {}). Exiting.\n'.format(bed_dem_file, num_rows, num_cols,\
    pixel_cell_map_file, num_rows_dem, num_cols_dem)

  # Read in the provided Bed Digital Elevation Map (BDEM) file to 2D bed_dem array
  logging.info('Loading Bed Digital Elevation Map (BDEM) from %s', bed_dem_file)
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
  # Read in the provided Surface Digital Elevation Map (SDEM) file to 2D 
  # surf_dem array
  logging.info('Loading Surface Digital Elevation Map (SDEM) from %s',\
    surf_dem_in_file)
  current_surf_dem = np.loadtxt(surf_dem_in_file, skiprows=5)

  # Check if Bed DEM has any points that are higher than the Surface DEM
  # in the same location.  If so, set these Bed DEM points to equal the
  # Surface DEM values, thus avoiding producing negative values when the
  # two are subtracted during glacier mask update. This reconciliation is
  # necessary because the two DEMs come from different sources, and could
  # have some overlapping elevation points.
  dem_diffs = current_surf_dem - bed_dem
  neg_val_inds = np.where(dem_diffs < 0)
  num_neg_vals = len(neg_val_inds)
  if num_neg_vals > 0:
    bed_dem[neg_val_inds] = current_surf_dem[neg_val_inds]
    new_bed_dem_file = bed_dem_file[0:-4] + '_adjusted.gsa'
    logging.warning('The provided Bed DEM (%s) has %s elevation points \
higher than those in the provided Surface DEM (%s), probably because they \
come from different data sources. The Bed DEM has been adjusted to equal \
the Surface DEM elevation at these points and written out to the file %s.',\
bed_dem_file, num_neg_vals, surf_dem_in_file, new_bed_dem_file)
    bed_dem_file = new_bed_dem_file
    write_grid_to_gsa_file(bed_dem, bed_dem_file, num_cols_dem, num_rows_dem,\
      0, 0, 0, 0)

  # Check agreement between elevation_map from VIC-grid-to-RGM-pixel file and
  # the loaded surface DEM. Mask out invalid values in elevation_map and get
  # valid indices to spatially align with subset of current_surf_dem that
  # falls within the simulation domain.
  masked_elevation_map = np.ma.masked_invalid(elevation_map)
  aligned_indices = masked_elevation_map.nonzero()
  assert (np.array_equal(masked_elevation_map[~masked_elevation_map.mask],\
    current_surf_dem[aligned_indices])), 'Values mismatch \
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
  logging.info('Loading initial Glacier Mask from %s', init_glacier_mask_file)
  glacier_mask = np.loadtxt(init_glacier_mask_file, skiprows=5)

  # Apply the initial glacier mask and modify the band and glacier area
  # fractions accordingly
  logging.debug('Applying initial area fraction update.')
  update_area_fracs(cells, cell_areas, cell_id_map, num_snow_bands,
    current_surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

  # Set the VIC output state filename prefix (to be written to STATENAME
  # in the global file)
  state_filename_prefix = temp_files_path + 'vic_hydrocon_state'

  # The RGM will always output a DEM file of the same name (if running RGM for
  # a single year at a time)
  rgm_surf_dem_out_file = temp_files_path + 's_out_00001.grd'

# (initialisation done)

#### Run the coupled VIC-RGM model for the time range specified in the VIC
  # global parameters file
  time_iterator = run_ranges(global_parms.startdate,
                 global_parms.enddate,
                 global_parms.glacier_accum_startdate)
  for start, end in time_iterator:

    # Write temporary VIC parameter files
    temp_snb = temp_files_path + 'snb_temp_' + start.isoformat() + '.txt'
    logging.debug('Writing temporary snow band parameter file %s', temp_snb)
    save_snb_parms(cells, temp_snb)
    temp_vpf = temp_files_path + 'vpf_temp_' + start.isoformat() + '.txt'
    logging.debug('Writing temporary vegetation parameter file %s', temp_vpf)
    save_veg_parms(cells, temp_vpf)
    temp_gpf = temp_files_path + 'gpf_temp_{}.txt'.format(start.isoformat())
    logging.debug('Writing temporary global parameter file %s', temp_gpf)
    global_parms.vegparam = temp_vpf
    global_parms.snow_band = '{} {}'.format(num_snow_bands, temp_snb)
    global_parms.startdate = start
    global_parms.enddate = end
    global_parms.statedate = end
    global_parms.statename = state_filename_prefix
    global_parms.write(temp_gpf)

    # Run VIC for a year, saving model state at the end
    print('\nRunning VIC from {} to {}'.format(start, end))
    logging.info('\nRunning VIC from %s to %s using global parameter file %s',\
      start, end, temp_gpf)
    try:
      subprocess.check_call([vic_full_path, "-g", temp_gpf], shell=False,\
        stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
      logging.error('Subprocess invocation of VIC failed with the following \
error: %s', e)
      sys.exit(0)

    # Open VIC NetCDF state file and load the most recent set of state
    # variable values for all grid cells being modeled
    state_file = state_filename_prefix + "_" + end.isoformat()
    logging.info('Reading saved VIC state file %s', state_file)
    # leave the state file open for modification later
    state_dataset = netCDF4.Dataset(state_file, 'r+')
    state_dataset.set_auto_mask(False)
    # these should never change within a run of the Hydro-Conductor
    Cell.Nlayers = state_dataset.state_nlayer
    Cell.Nnodes = state_dataset.state_nnode
    # drop unused fit error term from glacier mass balance polynomial
    Cell.NglacMassBalanceEqnTerms = state_dataset.state_nglac_mass_balance_eqn_terms - 1
    state = state_dataset.variables
    # read new states of all cells
    read_state(state, cells)
    # optionally leave the last VIC state file on disk 
    if not output_trace_files:
        os.remove(state_file)
    # pull Glacier Mass Balance polynomials out of cell states
    gmb_polys = {}
    cell_ids = []
    for cell_id in cells:
        cell_ids.append(cell_id)
        # leave off the "fit error" term at the end of the GMB polynomial
        # TODO: generalize this to handle variable length GMB polynomials
        gmb_polys[cell_id] = cells[cell_id].cell_state.variables\
          ['GLAC_MASS_BALANCE_EQN_TERMS'][0:3]

    # Translate mass balances using grid cell GMB polynomials and current
    # surface DEM into a 2D RGM mass balance grid (MBG)
    logging.debug('Converting glacier mass balance polynomials to 2D grid.')
    mass_balance_grid = mass_balances_to_rgm_grid(gmb_polys, cell_id_map,\
      current_surf_dem, num_rows_dem, num_cols_dem, cell_ids)
    # write Mass Balance Grid to ASCII file to direct the RGM to use as input
    mbg_file = temp_files_path + 'mass_balance_grid_' + end.isoformat()\
      + '.gsa'
    logging.debug('Writing RGM input glacier mass balance grid to file %s', mbg_file)
    write_grid_to_gsa_file(mass_balance_grid, mbg_file, num_cols_dem,\
      num_rows_dem, dem_xmin, dem_xmax, dem_ymin, dem_ymax)

    # Run RGM for one year, passing it the MBG, BDEM, SDEM
    logging.info('Running RGM for current year with parameter file %s, \
Bed DEM file %s, Surface DEM file %s, Mass Balance Grid file %s',\
      rgm_params_file, bed_dem_file, surf_dem_in_file, mbg_file)
    try:
      subprocess.check_call([rgm_full_path, "-p", rgm_params_file, "-b",\
        bed_dem_file, "-d", surf_dem_in_file, "-m", mbg_file, "-o",\
        temp_files_path, "-s", "0", "-e", "0" ], shell=False,\
        stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
      logging.error('Subprocess invocation of RGM failed with the following \
error: %s', e)
      sys.exit(0)

    # remove temporary files if not saving for offline inspection
    if not output_trace_files:
      os.remove(mbg_file)
      os.remove(rgm_surf_dem_out_file)

    # Read in new Surface DEM file from RGM output
    logging.debug('Reading Surface DEM file from RGM output %s',\
      rgm_surf_dem_out_file)
    current_surf_dem = np.loadtxt(rgm_surf_dem_out_file, skiprows=5)
    temp_surf_dem_file = temp_files_path + 'rgm_surf_dem_out_'\
      + start.isoformat() + '.gsa'
    os.rename(rgm_surf_dem_out_file, temp_surf_dem_file)
    # this will be fed back into RGM on next time step:
    surf_dem_in_file = temp_surf_dem_file

    # Update glacier mask
    logging.debug('Updating Glacier Mask')
    glacier_mask = update_glacier_mask(current_surf_dem, bed_dem,\
      num_rows_dem, num_cols_dem)
    if output_trace_files:
      glacier_mask_file = temp_files_path + 'glacier_mask_'\
        + start.isoformat() + '.gsa'
      write_grid_to_gsa_file(glacier_mask, glacier_mask_file, num_cols_dem,\
      num_rows_dem, dem_xmin, dem_xmax, dem_ymin, dem_ymax)
    
    # Update HRU and band area fractions and state for all VIC grid cells
    logging.debug('Updating VIC grid cell area fracs and states')
    update_area_fracs(cells, cell_areas, cell_id_map, num_snow_bands,\
      current_surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

    # Update the VIC state file with new state information
    logging.debug('Updating VIC state file')
    write_state(state, cells)
    state_dataset.close()

    # Set the new state file name VIC will have to read in on next iteration
    global_parms.init_state = state_filename_prefix + "_" + end.isoformat()

# Main program invocation.
if __name__ == '__main__':
  main()
