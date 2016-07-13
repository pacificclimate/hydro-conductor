"""io.py 

  This module contains utility functions related to file input/output for the
  vic_rgm_conductor.

"""

import numpy as np
import csv

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
    nx = int(headers['NCOLS'])
    ny = int(headers['NROWS'])
    # create an empty two dimensional array
    cell_id_map = np.empty((ny, nx))
    cell_id_map.fill(np.nan)
    _ = f.readline() # Consume the column headers
    for line in f:
      # Note that we are ignoring the median elevation data (5th column)
      _, i, j, _, _, cell_id = line.split()
      i, j = int(i), int(j)
      if cell_id != 'NA': #otherwise we leave it as np.NaN
        cell_id_map[i,j] = cell_id
      # Increment the pixel-granularity area within the grid cell
      if cell_id in cell_areas:
        cell_areas[cell_id] += 1
      else:
        cell_areas[cell_id] = 1
  vic_cell_mask = np.ma.masked_array(np.int32(cell_id_map))
  vic_cell_mask[np.where(np.isnan(cell_id_map))] = np.ma.masked

  return vic_cell_mask, cell_areas, nx, ny

def mass_balances_to_rgm_grid(gmb_polys, vic_cell_mask, surf_dem, bed_dem, \
  num_rows_dem, num_cols_dem):
  """ Translate mass balances from grid cell GMB polynomials to 2D RGM pixel \
    grid to use as one of the inputs to RGM
  """
  mass_balance_grid = np.ma.empty(vic_cell_mask.shape)
  mass_balance_grid[np.where(vic_cell_mask.mask)] = 0
  try:
    for row in range(num_rows_dem):
      for col in range(num_cols_dem):
        # only grab elevation for pixels that fall within a VIC cell
        if vic_cell_mask[row][col] is not np.ma.masked:
          cell_id = str(vic_cell_mask[row][col])
          # read most recent median elevation of this pixel
          median_elev = surf_dem[row][col]
          mass_balance_grid[row][col] = gmb_polys[cell_id][0] + median_elev \
            * (gmb_polys[cell_id][1] + median_elev * gmb_polys[cell_id][2])
        else:
          surf_dem[row][col] = bed_dem[row][col]
  except Exception as e:
    print('mass_balances_to_rgm_grid: Exception while processing pixel at \
row {} column {}: \n{}'.format(row, col, e))
    logging.error('mass_balances_to_rgm_grid: Exception while processing pixel \
at row %s column %s: \n %s', row, col, e)
    sys.exit(0)
  return mass_balance_grid

def read_gsa_headers(dem_file):
  """ Opens and reads the header metadata from a GSA Digital Elevation Map
    file, verifies agreement with the VIC-RGM mapping file metadata, and
    returns the x and y extents metadata
  """
  with open(dem_file, 'r') as f:
    # First line
    first_line = f.readline()
    assert first_line.startswith('DSAA'), 'read_gsa_headers({}): DSAA header \
      on first line of DEM file was not found or is malformed.  DEM file does \
      not conform to ASCII grid format.'.format(dem_file)
    # Second line
    num_cols, num_rows = f.readline().split()
    xmin, xmax = f.readline().split()
    ymin, ymax = f.readline().split()
    out_1 = [float(n) for n in (xmin, xmax, ymin, ymax)]
    out_2 = [int(x) for x in (num_rows, num_cols)]
  return out_1 + out_2

def write_grid_to_gsa_file(grid, outfilename, num_cols_dem, num_rows_dem,\
    dem_xmin, dem_xmax, dem_ymin, dem_ymax):
  """ Writes a 2D grid to ASCII file in the input format expected by the RGM \
  for DEM and mass balance grids """
  zmin = np.min(grid)
  zmax = np.max(grid)
  header_rows = [['DSAA'], [num_cols_dem, num_rows_dem], [dem_xmin, dem_xmax],\
    [dem_ymin, dem_ymax], [zmin, zmax]]
  with open(outfilename, 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ')
    for header_row in header_rows:
      writer.writerow(header_row)
    for row in grid:
      writer.writerow(row)

def read_state(state_in, cells):
  """Reads the most recent state variables from the VIC state file produced by
    the most recent VIC run and updates the CellState and HruState object
    members of each cell.
  """
  num_lons = len(state_in['lon'])
  def get_2D_cell_indices(count):
    """Returns the 2D lat/lon indices of the cell for accessing it from the
      state file
    """
    return count // num_lons, count % num_lons

  cell_idx = 0
  for cell_id, cell in cells.items():
    cell_lat_idx, cell_lon_idx = get_2D_cell_indices(cell_idx)
    cell_hru_idx = 0
    # read all cell state variables
    for variable in cell.cell_state.variables:
      if variable == 'lat':
        cell.cell_state.variables[variable] = state_in[variable][cell_lat_idx]
      elif variable == 'lon':
        cell.cell_state.variables[variable] = state_in[variable][cell_lon_idx]
      else:
        cell.cell_state.variables[variable] = \
          state_in[variable][cell_lat_idx][cell_lon_idx]
    for band in cell.bands:
      # HRUs are sorted by ascending veg_type_num in VIC state file
      for hru_veg_type in band.hru_keys_sorted:
        # read all HRU state variables with dimensions (lat, lon, hru)
        for variable in band.hrus[hru_veg_type].hru_state.variables:
          band.hrus[hru_veg_type].hru_state.variables[variable] = \
            state_in[variable][cell_lat_idx][cell_lon_idx][cell_hru_idx]
        cell_hru_idx += 1
    cell_idx += 1

def write_state(cells, old_dataset, new_dataset, new_state_date):
  """Takes the dataset from the last VIC state file, copies its static
    metadata and writes a new state file with static metadata, new
    dynamic metadata, and the new state variable values from the CellState and
    HruState object members of each Cell object in cells.
  """
  num_lons = len(old_dataset.variables['lon'])
  def get_2D_cell_indices(count):
    """Returns the 2D lat/lon indices of the cell for accessing it from the
      state file
    """
    return count // num_lons, count % num_lons

  new_dataset.state_year = np.int32(new_state_date.year)
  new_dataset.state_month = np.int32(new_state_date.month)
  new_dataset.state_day = np.int32(new_state_date.day)

  # Copy static attributes
  for attr in old_dataset.ncattrs():
    if attr not in ['state_year', 'state_month', 'state_day']:
      new_dataset.setncattr(attr, old_dataset.getncattr(attr))
  # Copy dimensions
  for d_name, dim in old_dataset.dimensions.items():
    if d_name == 'hru': # need to update HRU dimension because it can change
      max_num_hrus = max(sum([band.num_hrus for band in cell.bands]) for cell_id, cell in cells.items())
      new_dataset.createDimension(d_name, max_num_hrus)
    else:
      new_dataset.createDimension(d_name, len(dim) if not dim.isunlimited() else None)
  # Copy variables
  for v_name, var in old_dataset.variables.items():
    new_var = new_dataset.createVariable(v_name, var.datatype, var.dimensions)
    # Copy variable attributes
    new_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})

  state = new_dataset.variables
  cell_idx = 0
  for cell_id, cell in cells.items():
    cell_lat_idx, cell_lon_idx = get_2D_cell_indices(cell_idx)
    cell_hru_idx = 0
    # write all cell state variables
    for variable in cell.cell_state.variables:
      if variable == 'lat':
        state[variable][cell_lat_idx] = cell.cell_state.variables[variable]
      elif variable == 'lon':
        state[variable][cell_lon_idx] = cell.cell_state.variables[variable]
      else:
        state[variable][cell_lat_idx, cell_lon_idx] = \
          cell.cell_state.variables[variable]
    for band in cell.bands:
      # HRUs are sorted by ascending veg_type_num in VIC state file
      for hru_veg_type in band.hru_keys_sorted:
        # write all HRU state variables with dimensions (lat, lon, hru)
        for variable in band.hrus[hru_veg_type].hru_state.variables:
          state[variable][cell_lat_idx, cell_lon_idx, cell_hru_idx] = \
            band.hrus[hru_veg_type].hru_state.variables[variable]
        cell_hru_idx += 1
    cell_idx += 1
