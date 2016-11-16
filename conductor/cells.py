"""cells.py

  This module captures the conceptual components of the VIC grid cell,
  broken into classes by elevation band and Hydrological Response Unit (HRU).
  It also provides functions that operate on these class objects, nested
  within the OrderedDict of VIC cells. 
"""

from collections import OrderedDict
from copy import deepcopy
from math import ceil
import numpy as np
import itertools
import logging

from conductor.conductor_params import MAX_SURFACE_SWE, CH_ICE, \
  ZERO_AREA_FRAC_TOL, GLACIER_THICKNESS_THRESHOLD

# This is necessary pre-Python 3.5, after which point use math.isclose()
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

class Cell(object):
  """Class capturing VIC cells
  """
  # Dimensions Nlayers, Nnodes, dist, and NglacMassBalanceEqnTerms 
  # are defined on a *per-run* basis
  Nlayers = 3
  Nnodes = 3
  dist = 1
  NglacMassBalanceEqnTerms = 3

  def __init__(self, bands):
    self.bands = bands
    self.cell_state = CellState()

  def __eq__(self, other):
    return (self.__class__ == other.__class__ and self.__dict__ == other.__dict__)

  def update_cell_state(self):
    self.cell_state.variables['VEG_TYPE_NUM'] = sum([i.num_hrus for i in self.bands])
    # Nothing is to be done with the other CellState.variables; their values
    # that were read in by read_state() do not get modified.

  @property
  def num_bands(self):
    return len(self.bands)

class CellState(object):
  """Class capturing the set of VIC cell state and metadata variables that can
    change in a yearly VIC run.
  """
  def __init__(self):
    self.variables = {
      'lat':0.0,
      'lon':0.0,
      'GRID_CELL':0,
      'NUM_BANDS':0,
      'SOIL_DZ_NODE': [0]*Cell.Nnodes,
      'SOIL_ZSUM_NODE': [0]*Cell.Nnodes,
      'VEG_TYPE_NUM': 0,
      'GLAC_MASS_BALANCE_EQN_TERMS': [0]*Cell.NglacMassBalanceEqnTerms
    }

  def __repr__(self):
    return '{} (\n  '.format(self.__class__.__name__) + ' \n  '\
      .join([': '.join([key, str(value)]) \
        for key, value in self.variables.items()]) + '\n)'

  def __eq__(self, other):
    return (self.__class__ == other.__class__ and self.__dict__ == other.__dict__)

class Band(object):
  """Class capturing VIC cell parameters at the elevation band level
    (aka "snow band")
  """
  # glacier_id and open_ground_id, glacier_root_zone_parms, 
  # open_ground_root_zone_parms, and band_size are defined on a *per-run*
  # basis.
  glacier_id = 22
  glacier_root_zone_parms = [0.10, 1.00, 0.10, 0.00, 0.10, 0.00]
  open_ground_id = 19
  open_ground_root_zone_parms = [0.10, 1.00, 0.10, 0.00, 0.10, 0.00]
  band_size = 100

  def __init__(self, median_elev, hrus=None):
    self.median_elev = median_elev
    if hrus is None:
      hrus = {}
    self.hrus = hrus

  def __eq__(self, other):
    return (self.__class__ == other.__class__ and self.__dict__ == other.__dict__)

  @property
  def hru_keys_sorted(self):
    return sorted([int(k) for k in self.hrus.keys()])

  @property
  def lower_bound(self):
    return self.median_elev - self.median_elev % self.band_size

  @property
  def upper_bound(self):
    return self.lower_bound + self.band_size

  @property
  def num_hrus(self):
    return len(self.hrus)

  @property
  def area_frac(self): 
    """The Band area fraction, equal to the sum of HRU area fractions within
      this band, which should be equal to the total area fraction of the band
      as represented in the Snow Band Parameters file (initial conditions)
    """
    if self.hrus:
      return sum([hru.area_frac for hru in self.hrus.values()])
    else:
      return 0

  @property
  def area_frac_glacier(self):
    if self.glacier_id in self.hrus:
      return self.hrus[self.glacier_id].area_frac
    else:
      return 0

  @property
  def area_frac_non_glacier(self):
    return self.area_frac - self.area_frac_glacier

  @property
  def area_frac_open_ground(self):
    return sum([hru.area_frac for veg_type, hru in self.hrus.items()\
      if veg_type == self.open_ground_id])

  def create_hru(self, band_id, veg_type, area_frac):
    """Creates a new HRU of provided veg_type and area_frac
    """
    # Append new hru to existing dict of HRUs for this band
    if veg_type == self.glacier_id:
      self.hrus[veg_type] = HydroResponseUnit(area_frac,\
        self.glacier_root_zone_parms, band_id, veg_type)
    elif veg_type == self.open_ground_id:
      self.hrus[veg_type] = HydroResponseUnit(area_frac,\
        self.open_ground_root_zone_parms, band_id, veg_type)

  def delete_hru(self, veg_type):
    """Deletes an HRU of veg_type within the Band
    """
    del self.hrus[veg_type]

  def __del__(self):
    """ delete Band """
    pass

  def __repr__(self):
    return '{}({}, {} {})'.format(self.__class__.__name__, self.area_frac,\
      self.median_elev, repr(self.hrus))
  def __str__(self):
    return '{}({:.2f}% @{} meters with {} HRUs)'.\
      format(self.__class__.__name__, self.area_frac*100, self.median_elev,\
      len(self.hrus))

class HydroResponseUnit(object):
  """Class capturing vegetation parameters at the single vegetation
    tile (HRU) level (of which there can be many per band).
  """
  def __init__(self, area_frac, root_zone_parms, band_id, veg_type):
    self.area_frac = area_frac
    self.root_zone_parms = root_zone_parms
    self.hru_state = HruState(band_id, veg_type)
  def __repr__(self):
    return '{}({}, {})'.format(self.__class__.__name__,
                   self.area_frac, self.root_zone_parms)
  def __str__(self):
    return 'HRU({:.2f}%, {})'.format(self.area_frac*100, self.root_zone_parms)

  def __eq__(self, other):
    return (self.__class__==other.__class__ and self.__dict__==other.__dict__)

  def __ne__(self, other):
    return not self.__eq__(other)

class HruState(object):
  """Class capturing the set of VIC HRU state variables.
  """
  def __init__(self, band_id, veg_type):
    # variables is an OrderedDict because there is temporal dependence in the
    # state update among some of them when update_hru_state() is called
    self.variables = OrderedDict([
      # HRU state variables with dimensions (lat, lon, hru)
      ('HRU_BAND_INDEX', band_id),
      ('HRU_VEG_INDEX', veg_type),
      # These two have dimensions (lat, lon, hru, dist, Nlayers)
      ('LAYER_ICE_CONTENT', Cell.dist * [[0]*Cell.Nlayers]),
      ('LAYER_MOIST', Cell.dist * [[0]*Cell.Nlayers]),
      # HRU_VEG_VAR_WDEW has dimensions (lat, lon, hru, dist)
      ('HRU_VEG_VAR_WDEW' , [0]),
      # HRU state variables with dimensions (lat, lon, hru)
      ('SNOW_SWQ', 0),
      ('SNOW_DEPTH', 0),
      ('SNOW_DENSITY', 0),
      ('SNOW_CANOPY', 0),
      ('SNOW_PACK_WATER', 0),
      ('SNOW_SURF_WATER', 0),
      ('GLAC_WATER_STORAGE', 0),
      ('GLAC_CUM_MASS_BALANCE', 0),
      # HRU state variables with dimensions (lat, lon, hru, Nnodes)
      ('ENERGY_T', [0]*Cell.Nnodes),
      ('ENERGY_T_FBCOUNT', [0]*Cell.Nnodes),
      # HRU state variables with dimensions (lat, lon, hru)
      ('ENERGY_TFOLIAGE', 0),
      ('GLAC_SURF_TEMP', 0),
      ('SNOW_SURF_TEMP', 0),
      ('SNOW_COLD_CONTENT', 0),
      ('SNOW_PACK_TEMP', 0),
      ('SNOW_ALBEDO', 0),
      ('SNOW_LAST_SNOW', 0),
      ('SNOW_MELTING', 0),
      ('ENERGY_TFOLIAGE_FBCOUNT', 0),
      ('ENERGY_TCANOPY_FBCOUNT', 0),
      ('ENERGY_TSURF_FBCOUNT', 0),
      ('GLAC_SURF_TEMP_FBCOUNT', 0),
      ('SNOW_SURF_TEMP_FBCOUNT', 0),
      # remaining state variables from the "miscellaneous" list (lat, lon, hru)
      ('GLAC_SURF_TEMP_FBFLAG', 0),
      ('GLAC_VAPOR_FLUX', 0),
      ('SNOW_CANOPY_ALBEDO', 0),
      ('SNOW_SURFACE_FLUX', 0),
      ('SNOW_SURF_TEMP_FBFLAG', 0),
      ('SNOW_TMP_INT_STORAGE', 0),
      ('SNOW_VAPOR_FLUX', 0)
    ])
    
  def __repr__(self):
    return '{} (\n  '.format(self.__class__.__name__) + ' \n  '\
      .join([': '.join([key, str(value)]) \
        for key, value in self.variables.items()]) + '\n)'

  def __eq__(self, other):
    return (self.__class__ == other.__class__ and self.__dict__ == other.__dict__)

  def __ne__(self, other):
    return not self.__eq__(other)

# Following are the state variables split into sets according to their update
# method specification, as detailed in the VIC State Updating Spec 3.0.
# Note that the order of variables within the follow lists matters in some
# cases, where the updated value of variable is derived from the updated
# value of another.
spec_1_vars = ['SOIL_DZ_NODE', 'SOIL_ZSUM_NODE',\
  'VEG_TYPE_NUM']

spec_2_vars = ['LAYER_ICE_CONTENT', 'LAYER_MOIST',\
  'HRU_VEG_VAR_WDEW', 'SNOW_CANOPY', 'SNOW_DEPTH',\
  'SNOW_PACK_WATER', 'SNOW_SURF_WATER', 'SNOW_SWQ',\
  'SNOW_PACK_TEMP', 'SNOW_SURF_TEMP']

spec_3_vars = ['SNOW_DENSITY']

spec_4_vars = ['GLAC_WATER_STORAGE']

spec_5_vars = ['GLAC_CUM_MASS_BALANCE']

spec_6_vars = ['SNOW_COLD_CONTENT']

spec_7_vars = ['SNOW_ALBEDO', 'SNOW_LAST_SNOW', 'SNOW_MELTING']

spec_9_vars = ['ENERGY_T', 'ENERGY_TFOLIAGE', 'ENERGY_T_FBCOUNT',\
  'ENERGY_TCANOPY_FBCOUNT', 'ENERGY_TFOLIAGE_FBCOUNT','ENERGY_TSURF_FBCOUNT',\
  'GLAC_SURF_TEMP', 'GLAC_SURF_TEMP_FBCOUNT', 'SNOW_SURF_TEMP_FBCOUNT',\
  # miscellaneous vars:
  'GLAC_SURF_TEMP_FBFLAG', 'GLAC_VAPOR_FLUX',\
  'SNOW_CANOPY_ALBEDO', 'SNOW_SURFACE_FLUX', 'SNOW_SURF_TEMP_FBFLAG',\
  'SNOW_TMP_INT_STORAGE', 'SNOW_VAPOR_FLUX']

def apply_custom_root_zone_parms(hru_cell_dict, glacier_root_zone_parms,\
  open_ground_root_zone_parms):
  """Utility function to apply user-supplied custom root zone parameters
    to glacier and/or open ground HRUs at initialization time.
  """
  for cell in hru_cell_dict.keys():
    for key in hru_cell_dict[cell]:
      if glacier_root_zone_parms and (key[1] == global_parms.glacier_id):
        hru_cell_dict[cell][key].root_zone_parms = glacier_root_zone_parms
      if open_ground_root_zone_parms and (key[1] == global_parms.open_ground_id):
        hru_cell_dict[cell][key].root_zone_parms = open_ground_root_zone_parms

def merge_cell_input(hru_cell_dict, elevation_cell_dict):
  """Utility function to merge the dict of HRUs loaded via
    vegparams.load_veg_parms() with the list of Bands loaded at start-up via
    snbparams.load_snb_parms() into one unified structure capturing all VIC
    cells' initial properties (but not state).
  """ 
  missing_keys = hru_cell_dict.keys() ^ elevation_cell_dict.keys()
  if missing_keys:
    raise Exception("One or more cell IDs were found in one input file,"
        "but not the other. IDs: {}".format(missing_keys))
  # initialize new cell container
  cells = OrderedDict()
  for cell_id in elevation_cell_dict:
    cells[cell_id] = Cell(deepcopy(elevation_cell_dict[cell_id]))
  # FIXME: this is a little awkward
  for cell_id, hru_dict in hru_cell_dict.items():
    band_ids = { band_id for band_id, _ in hru_dict.keys() }
    band_dict = { band_id: {} for band_id in band_ids }
    for (band_id, veg_type), hru in hru_dict.items():
      band_dict[band_id][veg_type] = hru
    for band_id, hru_dict in band_dict.items():
      cells[cell_id].bands[band_id].hrus = hru_dict
    cells[cell_id].update_cell_state()
  return cells

def update_glacier_mask(surf_dem, bed_dem, num_rows_dem, num_cols_dem):
  """ Takes output Surface DEM from RGM and uses element-wise differencing 
    with the Bed DEM to form an updated glacier mask 
  """
  diffs = surf_dem - bed_dem
  if np.any(diffs < 0):
    raise Exception(
      'update_glacier_mask: Error: Subtraction of Bed DEM from the output \
      Surface DEM of RGM produced one or more negative values.'
    )

  glacier_mask = np.zeros((num_rows_dem, num_cols_dem))
  glacier_mask[diffs > GLACIER_THICKNESS_THRESHOLD] = 1
  return glacier_mask

def update_area_fracs(cells, cell_areas, vic_cell_mask, num_snow_bands,\
  surf_dem, num_rows_dem, num_cols_dem, glacier_mask, update_state=True):
  """Applies the updated RGM DEM and glacier mask and calculates and updates
    all HRU area fractions for all elevation bands within the VIC cells.
    Determines the HRU state update case based upon changes in HRU area
    fractions since the last time step, as per Algorithm Specification --
    VIC State Updating (Version 3.0.0), and calls update_hru_state().
  """
  def reverse_enumerate(iterable):
    """
    Enumerate over an iterable in reverse order while retaining proper indexes
    """
    return itertools.zip_longest(reversed(range(len(iterable))), reversed(iterable))

  # Create temporary band area and glacier area bins.
  # Counting pixels is a proxy for area within each band:
  band_areas = {}
  # Counting pixels landing within the glacier mask is a proxy for glacier area:
  glacier_areas = {}

  # Scan through the RGM output DEM and drop pixels into bins according to
  # cell and elevation
  for cell_id, cell in cells.items():
    logging.debug('Binning RGM pixels and updating area fractions for cell %s', cell_id)
    band_areas[cell_id] = [0] * num_snow_bands
    glacier_areas[cell_id] = [0] * num_snow_bands

    # Select the portion of the DEM that pertains to this cell from which we
    # will do binning
    masked_dem = np.ma.masked_array(surf_dem)
    masked_dem[np.where(vic_cell_mask != float(cell_id))] = np.ma.masked
    # Create a regular 'flat' np.array of the DEM subset that is within the
    # VIC cell
    flat_dem = masked_dem[~masked_dem.mask]

    # Check if any pixels fall outside of valid range of bands
    if len(np.where(flat_dem < cell.bands[0].lower_bound)[0]) > 0:
      raise Exception(
        'One or more RGM output DEM pixels lies below the bounds of the lowest '\
        'defined elevation band (< {}m) as defined by the Snow Band Parameter File '\
        'for cell {}. You may need to add or shift the zero padding to accommodate this.'\
        .format(cell.bands[0].lower_bound, cell_id))
    if len(np.where(flat_dem >= cell.bands[-1].upper_bound)[0]) > 0:
      raise Exception(
        'One or more RGM output DEM pixels lies above the bounds of the highest '\
        'defined elevation band (>= {}m) as defined by the Snow Band Parameter File '\
        'for cell {}. You may need to add or shift the zero padding to accommodate this.'\
        .format(cell.bands[-1].upper_bound, cell_id))

    # Identify the bounds of the band area and glacier area bins for this cell
    # and update the band median elevations 
    band_bin_bounds = []
    for band in cell.bands:
      band_bin_bounds.append(band.lower_bound)
      band_pixels = flat_dem[np.where((flat_dem >= band.lower_bound) &\
        (flat_dem < band.upper_bound))]
      if not band_pixels.size: # if there are no pixels in this band
        band.median_elev = band.lower_bound
      else:
        band.median_elev = np.median(band_pixels)

    # Select portion of glacier_mask pertaining to this cell
    cell_glacier_mask = np.ma.masked_array(glacier_mask)
    cell_glacier_mask[np.where(vic_cell_mask != float(cell_id))] = np.ma.masked
    masked_dem.mask = False
    masked_dem[np.where(cell_glacier_mask != 1)] = np.ma.masked
    # Create a regular 'flat' np.array of the DEM subset that is unmasked glacier
    flat_glacier_dem = masked_dem[~masked_dem.mask]

    # Do binning into band_areas[cell][band] and glacier_areas[cell][band] 
    # using data in flat_dem and flat_glacier_dem
    inds = np.digitize(flat_dem, band_bin_bounds)
    for idx, count in enumerate(np.bincount(inds-1)):
      band_areas[cell_id][idx] = count

    inds = np.digitize(flat_glacier_dem, band_bin_bounds)
    for idx, count in enumerate(np.bincount(inds-1)):
      glacier_areas[cell_id][idx] = count

    # Initialise temporary band-level area fractions used
    # in computations below
    new_band_area_frac = [0] * cell.num_bands
    new_glacier_area_frac = [0] * cell.num_bands
    new_non_glacier_area_frac = [0] * cell.num_bands
    new_open_ground_area_frac = [0] * cell.num_bands
    veg_scaling_divisor = [0] * cell.num_bands
    delta_area_vegetated = [0] * cell.num_bands
    delta_area_hru = {}
    new_hru_area_frac = {}

    ### Update all Band area fractions for this cell
    # This must be done for all bands before the HRU area
    # fractions and state update steps
    for band_id, band in reverse_enumerate(cell.bands):
      # Update total area fraction for this Band
      new_band_area_frac[band_id] = band_areas[cell_id][band_id]/cell_areas[cell_id]
      # Update glacier HRU area fraction for this Band
      new_glacier_area_frac[band_id] = glacier_areas[cell_id][band_id]/cell_areas[cell_id]
      if isclose(new_glacier_area_frac[band_id], 0, abs_tol=ZERO_AREA_FRAC_TOL):
        new_glacier_area_frac[band_id] = 0

      # We need to update HRU area fractions if either 
      # the glacier HRU area fraction has changed for this band, or
      # the band's total area fraction has changed
      if (new_glacier_area_frac[band_id] != band.area_frac_glacier)\
        or (new_band_area_frac[band_id] != band.area_frac):
        # Calculate new non-glacier area fraction for this band
        new_non_glacier_area_frac[band_id] = new_band_area_frac[band_id] - new_glacier_area_frac[band_id]
        # Calculate new_residual area fraction
        new_residual_area_frac = new_non_glacier_area_frac[band_id]\
          - band.area_frac_non_glacier
        # Calculate new open ground area fraction
        new_open_ground_area_frac[band_id] = np.max([0, (band.area_frac_open_ground\
          + new_residual_area_frac)])
        if isclose(new_open_ground_area_frac[band_id], 0, abs_tol=ZERO_AREA_FRAC_TOL):
          new_open_ground_area_frac[band_id] = 0
        # Use old proportions of vegetated areas for scaling their area
        # fractions in this iteration
        veg_scaling_divisor[band_id] = band.area_frac_non_glacier\
          - band.area_frac_open_ground
        # Calculate the change in sum of vegetated area fractions
        delta_area_vegetated[band_id] = np.min([0, (band.area_frac_open_ground\
          + new_residual_area_frac)])

        for veg_type, hru in band.hrus.items():
          if veg_type is not Band.glacier_id and veg_type is not Band.open_ground_id:
            # Calculate change in HRU area fraction & update
            delta_area_hru[str(veg_type)] = delta_area_vegetated[band_id] * \
              (band.hrus[veg_type].area_frac / veg_scaling_divisor[band_id])
            if str(band_id) not in new_hru_area_frac:
              new_hru_area_frac[str(band_id)] = {}
            new_hru_area_frac[str(band_id)][str(veg_type)] = \
              band.hrus[veg_type].area_frac + delta_area_hru[str(veg_type)]
            if isclose(new_hru_area_frac[str(band_id)][str(veg_type)], 0, abs_tol=ZERO_AREA_FRAC_TOL):
              new_hru_area_frac[str(band_id)][str(veg_type)] = 0

    if update_state:
      # Update all HRU states for each band, then apply new HRU area
      # fractions calculated above.
      for band_id, band in reverse_enumerate(cell.bands):
        update_band_state(cell, band, band_id, new_band_area_frac,
          new_glacier_area_frac, new_open_ground_area_frac, new_hru_area_frac,
          delta_area_hru)
      # Update cell-level state variables
      cell.update_cell_state()

        # TODO: replace test below with one that checks that all band HRU area fracs add to 1
        # Sanity check that glacier + non-glacier area fractions add up to Band's total area fraction, within tolerance
        # sum_test = new_glacier_area_frac[band_id] + new_non_glacier_area_frac[band_id]
        # if sum_test != new_band_area_frac[band_id]:
        #   raise Exception(
        #       'cell {}, band {}: glacier area fraction {} + '
        #       'non-glacier area fraction {} = {} is not equal to '
        #       'the band area fraction of {}'
        #       .format(cell_id, band_id, new_glacier_area_frac,
        #           new_non_glacier_area_frac[band_id], sum_test,
        #           new_band_area_frac)
        #   )

def update_band_state(cell, band, band_id, new_band_area_frac,
            new_glacier_area_frac, new_open_ground_area_frac, new_hru_area_frac,
            delta_area_hru):
  """
  Updates all HRU states for each band, then applies newly
  calculated HRU area fractions provided.
  """
  # HRU area fractions from the previous time step are held in the
  # Band / HRU members of a given Cell instance
  # (e.g. band.area_frac_glacier or band.hrus[veg_type].area_frac)
  # and are compared against their newly calculated values for the
  # current time step (e.g. new_glacier_area_frac[band_id]) to determine
  # whether HRUs should be added/deleted, and which state update case
  # should be carried out by update_hru_state().
  # The algorithm requires that we update glacier and open ground HRU
  # states first, followed by all remaining HRU vegetation types.

  previous_glacier_area_fracs = [0] * cell.num_bands

  hrus_to_be_deleted = []  # HRUs marked for deletion at end of band iteration

  # We need to update area fractions and update state in a band if
  # the band's total area fraction has changed, or the band's glacier
  # HRU area fraction has changed.
  logging.debug('Starting state update for band %s...', band_id)
  if (new_glacier_area_frac[band_id] != band.area_frac_glacier) \
          or (new_band_area_frac[band_id] != band.area_frac):

    ## Glacier HRU update:
    logging.debug('GLACIER HRU %s update phase...', Band.glacier_id)
    # NOTE: we never mark glacier HRUs for deletion, we just set their
    # area fractions to zero. This is VIC's 'shadow glacier' requirement.
    if new_glacier_area_frac[band_id] > 0 and band.area_frac_glacier == 0:
      # CASE 1: A new glacier HRU has appeared in this band. Create a new
      # glacier HRU (HRU state defaults are automatically set upon HRU
      # creation, so technically there's no need to call update_hru_state())
      logging.debug('State update CASE 1 identified. New glacier appeared. '
        'Creating GLACIER HRU with area fraction %s',
        new_glacier_area_frac[band_id])
      band.create_hru(band_id, Band.glacier_id, new_glacier_area_frac[band_id])
      new_area_fracs = {}
      update_hru_state(None, None, '1', **new_area_fracs)
    elif new_glacier_area_frac[band_id] == band.area_frac_glacier:
      # CASE 2: Trivial case. There was no change in glacier HRU area,
      # so we don't call update_hru_state())
      logging.debug('State update CASE 2 identified. No change in GLACIER '
        'HRU area. No state update required.')
    elif new_glacier_area_frac[band_id] > 0 and band.area_frac_glacier > 0 \
            and new_glacier_area_frac[band_id] != band.area_frac_glacier:
      # CASE 3: Glacier HRU exists in previous and current time steps, but
      # its area fraction has changed.
      logging.debug('State update CASE 3 identified. GLACIER HRU area '
        'fraction has changed from %s to %s.', band.area_frac_glacier,
        new_glacier_area_frac[band_id])
      new_area_fracs = {
        'new_hru_area_frac': new_glacier_area_frac[band_id]
      }
      update_hru_state(
        band.hrus[Band.glacier_id],
        band.hrus[Band.glacier_id],
        '3', **new_area_fracs)
    elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
            and new_band_area_frac[band_id] > 0 and band.area_frac > 0:
      # CASE 4a: Glacier HRU has disappeared, but the band remains
      # (implies open ground HRU is expanding). Add state to open ground
      # HRU (leave the zero area "shadow glacier" HRU in place for VIC).
      logging.debug('State update CASE 4a identified. Glacier has '
        'disappeared. Transferring state to the OPEN GROUND HRU in this band.')
      new_area_fracs = {
        'new_open_ground_area_frac': new_open_ground_area_frac[band_id]
      }
      # if there's not already an open ground HRU in this band, create one
      if Band.open_ground_id not in band.hrus:
        band.create_hru(band_id, Band.open_ground_id, new_open_ground_area_frac[band_id])
      update_hru_state(
        band.hrus[Band.glacier_id],
        band.hrus[Band.open_ground_id],
        '4a', **new_area_fracs)
    elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
            and new_band_area_frac[band_id] == 0 and band.area_frac != 0 \
            and (band_id - 1) >= 0 \
            and new_glacier_area_frac[band_id - 1] > 0:
      # CASE 5a: Glacier HRU and the band have disappeared.
      # An adjacent lower band exists and there's glacier in it
      # (with non-zero area fraction, i.e. not a shadow glacier). Add
      # state to the glacier HRU in the lower band.
      logging.debug('State update CASE 5a identified. Glacier and the band '
        'have disappeared. Transferring state to the GLACIER HRU in band %s '
        'below.', band_id - 1)
      if not Band.glacier_id in cell.bands[band_id - 1].hrus:
        # Due to the top-to-bottom iteration ordering of band processing,
        # we need to create a glacier HRU in the lower band if one doesn't
        # already exist
        logging.debug('(CASE 5a) New GLACIER has appeared in band %s. '
          'Creating GLACIER HRU with area fraction %s', band_id - 1,
          new_glacier_area_frac[band_id - 1])
        cell.bands[band_id - 1].create_hru(band_id - 1, Band.glacier_id, \
                                           new_glacier_area_frac[band_id - 1])
      new_area_fracs = {
        'new_glacier_area_frac': new_glacier_area_frac[band_id - 1]
      }
      update_hru_state(
        band.hrus[Band.glacier_id],
        cell.bands[band_id - 1].hrus[Band.glacier_id],
        '5a', **new_area_fracs)
    elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
            and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
            and (band_id - 1) >= 0 \
            and new_open_ground_area_frac[band_id - 1] > 0:
      # CASE 5b: Glacier HRU and the band have disappeared.
      # An adjacent lower band exists with a non-zero area fraction open
      # ground HRU in it. Add state to that lower band's open ground HRU.
      logging.debug('State update CASE 5b identified. Glacier and the band '
        'have disappeared. Transferring state to the OPEN GROUND HRU in '
        'band %s below.', band_id - 1)
      if not Band.open_ground_id in cell.bands[band_id - 1].hrus:
        # Due to the top-to-bottom iteration ordering of band processing,
        # we need to create an open ground HRU in the lower band if one
        # doesn't already exist
        logging.debug('(CASE 5b) New OPEN GROUND has appeared in band %s. '
          'Creating OPEN GROUND HRU with area fraction %s.', band_id - 1,
          new_open_ground_area_frac[band_id - 1])
        cell.bands[band_id - 1].create_hru(band_id - 1, Band.open_ground_id, \
                                           new_open_ground_area_frac[band_id - 1])
      new_area_fracs = {
        'new_open_ground_area_frac': new_open_ground_area_frac[band_id - 1]
      }
      update_hru_state(
        band.hrus[Band.glacier_id],
        cell.bands[band_id - 1].hrus[Band.open_ground_id],
        '5b', **new_area_fracs)
    elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
            and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
            and (band_id - 1) >= 0 and new_band_area_frac[band_id - 1] > 0:
      # CASE 5c: Glacier HRU and the band have disappeared due to glacier.
      # If an adjacent lower band exists (with non-zero area fraction), add
      # state to the vegetated HRU with non-zero area fraction and the
      # greatest vegetation type index in that band.
      valid_vegetated_hrus_below = list(filter(lambda x: x[1].area_frac > 0, \
                                               new_hru_area_frac[str(band_id - 1)].items()))
      max_veg_type_below = max(valid_vegetated_hrus_below)[0]
      logging.debug('State update CASE 5c identified. Glacier and the band '
        'have disappeared. Transferring state to the VEGETATED HRU %s in '
        'band %s below.', max_veg_type_below, band_id - 1)
      new_area_fracs = {
        'new_hru_area_frac': new_hru_area_frac[str(band_id - 1)][str(max_veg_type_below)]
      }
      update_hru_state(
        band.hrus[Band.glacier_id],
        cell.bands[band_id - 1].hrus[int(max_veg_type_below)],
        '5c', **new_area_fracs)
    elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
            and new_band_area_frac[band_id] == 0 and band.area_frac != 0 \
            and ((band_id - 1) < 0 or new_band_area_frac[band_id - 1] == 0):
      # CASE 5d: Both the glacier HRU and band have disappeared due
      # to glacier. No adjacent lower band exists (or not one has non-zero
      # area fraction) to transfer state to. In this special case, add
      # state to glacier in the band above.
      logging.debug('State update CASE 5d identified. Glacier and the band '
        'have disappeared, and no adjacent lower band exists. Transferring '
        'state to the GLACIER HRU in band %s above.', band_id + 1)
      new_area_fracs = {
        'new_glacier_area_frac': new_glacier_area_frac[band_id + 1]
      }
      update_hru_state(
        band.hrus[Band.glacier_id],
        cell.bands[band_id + 1].hrus[Band.glacier_id],
        '5d', **new_area_fracs)
    else:
      raise Exception(
        'Error: No state update case identified for cell {}, band {}, HRU {}.'
        .format(cell_id, band_id, Band.glacier_id)
      )

    # Save the glacier area fraction for this band for checking against
    # in the next (lower) band iteration (CASE 5d specifically)
    previous_glacier_area_fracs[band_id] = band.area_frac_glacier

    # Apply update to glacier HRU area fraction
    if Band.glacier_id in band.hrus:
      band.hrus[Band.glacier_id].area_frac = new_glacier_area_frac[band_id]

    ## Open ground HRU update:
    logging.debug('OPEN GROUND HRU %s update phase...', Band.open_ground_id)
    if new_open_ground_area_frac[band_id] > 0 and band.area_frac_open_ground == 0:
      # CASE 1: New open ground was exposed in this band, so
      # create a new open ground HRU (HRU state defaults are
      # automatically set upon HRU creation, so technically
      # there's no need to call update_hru_state())
      logging.debug('State update CASE 1 identified. New open ground '
        'exposed. Creating OPEN GROUND HRU with area fraction %s.',
        new_open_ground_area_frac[band_id])
      band.create_hru(band_id, Band.open_ground_id,
                      new_open_ground_area_frac[band_id])
      new_area_fracs = {}
      update_hru_state(None, None, '1', **new_area_fracs)
    elif new_open_ground_area_frac[band_id] == band.area_frac_open_ground:
      # CASE 2: Trivial case. There was no change in open ground HRU area,
      # so we don't call update_hru_state())
      logging.debug('State update CASE 2 identified. No change in OPEN '
        'GROUND HRU area. No state update required.')
    elif new_open_ground_area_frac[band_id] > 0 and band.area_frac_open_ground > 0 \
            and new_open_ground_area_frac[band_id] != band.area_frac_open_ground:
      # CASE 3: Open ground HRU exists in previous and current time
      # steps, but its area fraction has changed.
      logging.debug('State update CASE 3 identified. OPEN GROUND HRU '
        'area fraction has changed from %s to %s.',
        band.area_frac_open_ground, new_open_ground_area_frac[band_id])
      new_area_fracs = {
        'new_hru_area_frac': new_open_ground_area_frac[band_id]
      }
      update_hru_state(
        band.hrus[Band.open_ground_id],
        band.hrus[Band.open_ground_id],
        '3', **new_area_fracs)
    elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
            and new_band_area_frac[band_id] != 0 and band.area_frac > 0:
      # CASE 4b: Open ground has disappeared, but the band remains
      # (implies glacier HRU expanding). Add state to the glacier HRU.
      # Mark the open ground HRU for deletion.
      logging.debug('State update CASE 4b identified. Open ground has '
        'disappeared. Transferring state to the GLACIER HRU in this band.')
      hrus_to_be_deleted.append(Band.open_ground_id)
      new_area_fracs = {
        'new_glacier_area_frac': new_glacier_area_frac[band_id]
      }
      update_hru_state(
        band.hrus[Band.open_ground_id],
        band.hrus[Band.glacier_id],
        '4b', **new_area_fracs)
    elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
            and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
            and (band_id - 1) >= 0 \
            and new_glacier_area_frac[band_id - 1] > 0:
      # CASE 5a: Open ground HRU and the band have disappeared.
      # An adjacent lower band exists and there's glacier in it
      # (with non-zero area fraction, i.e. not a shadow glacier). Add
      # state to the glacier HRU in the lower band. Mark the open
      # ground HRU for deletion.
      logging.debug('State update CASE 5a identified. Open ground and the '
        'band have disappeared. Transferring state to the GLACIER HRU in band '
        '%s below.', band_id - 1)
      hrus_to_be_deleted.append(Band.open_ground_id)
      if not Band.glacier_id in cell.bands[band_id - 1].hrus:
        # Due to the top-to-bottom iteration ordering of band processing,
        # we need to create a glacier HRU in the lower band if one doesn't
        # already exist (only if new_glacier_area_frac[band_id - 1] > 0)
        logging.debug('(CASE 5a) New GLACIER has appeared in band %s. '
          'Creating GLACIER HRU with area_frac %s.',
          band_id - 1, new_glacier_area_frac[band_id - 1])
        cell.bands[band_id - 1].create_hru(band_id - 1, Band.glacier_id, \
                                           new_glacier_area_frac[band_id - 1])
      new_area_fracs = {
        'new_glacier_area_frac': new_glacier_area_frac[band_id - 1]
      }
      update_hru_state(
        band.hrus[Band.open_ground_id],
        cell.bands[band_id - 1].hrus[Band.glacier_id],
        '5a', **new_area_fracs)
    elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
            and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
            and (band_id - 1) >= 0 \
            and new_open_ground_area_frac[band_id - 1] > 0:
      # CASE 5b: Open ground HRU and the band have disappeared.
      # An adjacent lower band exists with a non-zero area fraction open
      # ground HRU in it. Add state to that lower band's open ground HRU.
      # Mark the open ground HRU for deletion.
      logging.debug('State update CASE 5b identified. Open ground and the '
        'band have disappeared. Transferring state to the OPEN GROUND HRU in '
        'band %s below.', band_id - 1)
      hrus_to_be_deleted.append(Band.open_ground_id)
      if not Band.open_ground_id in cell.bands[band_id - 1].hrus:
        # Due to the top-to-bottom iteration ordering of band processing,
        # we need to create an open ground HRU in the lower band if one
        # doesn't already exist
        logging.debug('(CASE 5b) New OPEN GROUND has appeared in band %s. '
          'Creating OPEN GROUND HRU with area fraction %s.', band_id - 1,
          new_open_ground_area_frac[band_id - 1])
        cell.bands[band_id - 1].create_hru(band_id - 1, Band.open_ground_id,
                                        new_open_ground_area_frac[band_id - 1])
      new_area_fracs = {
        'new_open_ground_area_frac': new_open_ground_area_frac[band_id - 1]
      }
      update_hru_state(
        band.hrus[Band.open_ground_id],
        cell.bands[band_id - 1].hrus[Band.open_ground_id],
        '5b', **new_area_fracs)
    elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
            and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
            and (band_id - 1) >= 0 and new_band_area_frac[band_id - 1] > 0:
      # CASE 5c: Open ground HRU and the band have disappeared due to
      # glacier. Mark the open ground HRU for deletion. If an adjacent
      # lower band exists (with non-zero area_frac), add state to the
      # vegetated HRU with non-zero area_frac and the greatest vegetation
      # type index in that band.
      valid_vegetated_hrus_below = list(filter(lambda x: x[1].area_frac > 0,
                                    new_hru_area_frac[str(band_id - 1)].items()))
      max_veg_type_below = max(valid_vegetated_hrus_below)[0]
      logging.debug('State update CASE 5c identified. Open ground and the '
        'band have disappeared. Transferring state to the VEGETATED HRU %s in '
        'band %s below.', max_veg_type_below, band_id - 1)
      hrus_to_be_deleted.append(Band.open_ground_id)
      new_area_fracs = {
        'new_hru_area_frac': new_hru_area_frac[str(band_id - 1)][str(max_veg_type_below)]
      }
      update_hru_state(
        band.hrus[Band.open_ground_id],
        cell.bands[band_id - 1].hrus[int(max_veg_type_below)],
        '5c', **new_area_fracs)
    elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
            and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
            and ((band_id - 1) < 0 or new_band_area_frac[band_id - 1] == 0):
      # CASE 5d: Both the open ground HRU and band have disappeared due
      # to glacier. Mark the HRU for deletion. No adjacent lower band
      # exists (or not one has non-zero area fraction) to transfer
      # state to. In this special case, add state to glacier in the
      # band above.
      logging.debug('State update CASE 5d identified. Open ground and the '
        'band have disappeared, and no adjacent lower band exists. Transferring '
        'state to the GLACIER HRU in band %s above.', band_id + 1)
      hrus_to_be_deleted.append(Band.open_ground_id)
      new_area_fracs = {
        'new_glacier_area_frac': new_glacier_area_frac[band_id + 1]
      }
      update_hru_state(
        band.hrus[Band.open_ground_id],
        cell.bands[band_id + 1].hrus[Band.glacier_id],
        '5d', **new_area_fracs)
    else:
      raise Exception(
        'Error: No state update case identified for cell {}, band {}, HRU {}.'
        .format(cell_id, band_id, Band.open_ground_id)
      )
    # Apply update to open ground HRU area fraction. HRU will get deleted later
    # if this is 0 (it will have already been added to hrus_to_be_deleted)
    if Band.open_ground_id in band.hrus:
      band.hrus[Band.open_ground_id].area_frac = new_open_ground_area_frac[band_id]

    ### Update area fractions and states for all remaining non-glacier and
    # non-open ground HRUs in this band
    for veg_type, hru in band.hrus.items():
      if veg_type is not Band.glacier_id and veg_type is not Band.open_ground_id:
        logging.debug('VEGETATED HRU %s update phase...', veg_type)
        if new_hru_area_frac[str(band_id)][str(veg_type)] < 0:
          raise Exception(
            'Error: cell {}, band {}: HRU {} has a negative area fraction ({}).'
            .format(cell_id, band_id, veg_type, \
                    new_hru_area_frac[str(band_id)][str(veg_type)])
          )
        if new_hru_area_frac[str(band_id)][str(veg_type)] == band.hrus[veg_type].area_frac:
          # CASE 2: Trivial case. There was no change in this HRU's area,
          # so we don't call update_hru_state())
          logging.debug('State update CASE 2 identified. No change in '
            'VEGETATED HRU %s area. No state update required.', veg_type)
        elif new_hru_area_frac[str(band_id)][str(veg_type)] > 0 \
                and delta_area_hru[str(veg_type)] != 0:
          # CASE 3: Vegetated HRU exists in previous and current time
          # steps, but its area fraction has changed.
          logging.debug('State update CASE 3 identified. HRU %s area '
            'fraction has changed from %s to %s.',
            veg_type, new_hru_area_frac[str(band_id)][str(veg_type)] - \
            delta_area_hru[str(veg_type)],
            new_hru_area_frac[str(band_id)][str(veg_type)])
          new_area_fracs = {
            'new_hru_area_frac': new_hru_area_frac[str(band_id)][str(veg_type)]
          }
          update_hru_state(
            hru,
            hru,
            '3', **new_area_fracs)
        elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 \
                and new_band_area_frac[band_id] > 0 and band.area_frac > 0:
          # CASE 4b: Vegetated HRU has disappeared, but the band remains
          # (implies glacier HRU expanding). Mark the HRU for deletion
          # (never to return -- only open ground can come back in its place)
          # and add state to the glacier HRU.
          logging.debug('State update CASE 4b identified. HRU %s has '
            'disappeared. Transferring state to the GLACIER HRU in this band.',
            veg_type)
          hrus_to_be_deleted.append(veg_type)
          new_area_fracs = {
            'new_glacier_area_frac': new_glacier_area_frac[band_id]
          }
          update_hru_state(
            hru,
            band.hrus[Band.glacier_id],
            '4b', **new_area_fracs)
        elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 and band.area_frac > 0 \
                and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
                and (band_id - 1) >= 0 \
                and new_glacier_area_frac[band_id - 1] > 0:
          # CASE 5a: Both the vegetated HRU and the band have disappeared.
          # An adjacent lower band exists and there's glacier in it
          # (with non-zero area fraction, i.e. not a shadow glacier). Add
          # state to the glacier HRU in the lower band. Mark the vegetated
          # HRU for deletion.
          logging.debug('State update CASE 5a identified. HRU %s and the '
            'band have disappeared. Transferring state to the GLACIER HRU in '
            'band %s below.', veg_type, band_id - 1)
          hrus_to_be_deleted.append(veg_type)
          if not Band.glacier_id in cell.bands[band_id - 1].hrus:
            # Due to the top-to-bottom iteration ordering of band processing,
            # we need to create a glacier HRU in the lower band if one doesn't
            # already exist (only if new_glacier_area_frac[band_id - 1] > 0)
            logging.debug('(CASE 5a) New GLACIER has appeared in band %s. '
              'Creating GLACIER HRU with area fraction %s.', band_id - 1,
              new_glacier_area_frac[band_id - 1])
            cell.bands[band_id - 1].create_hru(band_id - 1, Band.glacier_id, \
                                               new_glacier_area_frac[band_id - 1])
          new_area_fracs = {
            'new_glacier_area_frac': new_glacier_area_frac[band_id - 1]
          }
          update_hru_state(
            hru,
            cell.bands[band_id - 1].hrus[Band.glacier_id],
            '5a', **new_area_fracs)
        elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 and band.area_frac > 0 \
                and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
                and (band_id - 1) >= 0 \
                and new_open_ground_area_frac[band_id - 1] > 0:
          # CASE 5b: Both the vegetated HRU and the band have disappeared.
          # An adjacent lower band exists with a non-zero area fraction open
          # ground HRU in it. Add state to that lower band's open ground HRU.
          # Mark the vegetated HRU for deletion.
          logging.debug('State update CASE 5b identified. HRU %s and the band '
            'have disappeared. Transferring state to the OPEN GROUND HRU in '
            'band %s below.', veg_type, band_id - 1)
          hrus_to_be_deleted.append(veg_type)
          if not Band.open_ground_id in cell.bands[band_id - 1].hrus:
            # Due to the top-to-bottom iteration ordering of band processing,
            # we need to create an open ground HRU in the lower band if one
            # doesn't already exist
            logging.debug('(CASE 5b) New OPEN GROUND has appeared in band %s. '
              'Creating OPEN GROUND HRU with area fraction %s.',
              band_id - 1, new_open_ground_area_frac[band_id - 1])
            cell.bands[band_id - 1].create_hru(band_id - 1, Band.open_ground_id,
                                          new_open_ground_area_frac[band_id - 1])
          new_area_fracs = {
            'new_open_ground_area_frac': new_open_ground_area_frac[band_id - 1]
          }
          update_hru_state(
            hru,
            cell.bands[band_id - 1].hrus[Band.open_ground_id],
            '5b', **new_area_fracs)
        elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 \
                and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
                and (band_id - 1) >= 0 and new_band_area_frac[band_id - 1] > 0:
          # CASE 5c: Both the vegetated HRU and the band have disappeared due to
          # glacier. Mark the vegetated HRU for deletion. An adjacent
          # lower band exists (with non-zero area fraction), where we will
          # add state to the vegetated HRU with non-zero area fraction and
          # the greatest vegetation type index in that band.
          valid_vegetated_hrus_below = list(filter(lambda x: x[1].area_frac > 0,
                                                   new_hru_area_frac[str(band_id - 1)].items()))
          max_veg_type_below = max(valid_vegetated_hrus_below)[0]
          logging.debug('State update CASE 5c identified. HRU %s and the '
            'band have disappeared. Transferring state to the VEGETATED HRU %s '
            'in band %s below.', veg_type, max_veg_type_below, band_id - 1)
          hrus_to_be_deleted.append(veg_type)
          new_area_fracs = {
            'new_hru_area_frac': new_hru_area_frac[str(band_id - 1)][str(max_veg_type_below)]
          }
          update_hru_state(
            hru,
            cell.bands[band_id - 1].hrus[int(max_veg_type_below)],
            '5c', **new_area_fracs)
        elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 \
                and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
                and ((band_id - 1) < 0 or new_band_area_frac[band_id - 1] == 0):
          # CASE 5d: Both the vegetated HRU and band have disappeared due
          # to glacier. Mark the vegetated HRU for deletion. No adjacent
          # lower band exists (or not one has non-zero area fraction) to
          # transfer state to. In this special case, add state to glacier
          # in the band above.
          logging.debug('State update CASE 5d identified. HRU %s and the '
            'band have disappeared, and no adjacent lower band exists. '
            'Transferring state to the GLACIER HRU in band %s above.',
            veg_type, band_id + 1)
          hrus_to_be_deleted.append(veg_type)
          new_area_fracs = {
            'new_glacier_area_frac': new_glacier_area_frac[band_id + 1]
          }
          update_hru_state(
            hru,
            cell.bands[band_id + 1].hrus[Band.glacier_id],
            '5d', **new_area_fracs)
        else:
          raise Exception(
            'Error: No state update case identified for cell {}, band {}, HRU {}.'
            .format(cell_id, band_id, veg_type)
          )
        # Apply update to the HRU's area fraction. HRU will get deleted later
        # if this is 0 (it will have already been added to hrus_to_be_deleted)
        band.hrus[veg_type].area_frac = new_hru_area_frac[str(band_id)][str(veg_type)]

    # Remove HRUs marked for deletion
    for hru in hrus_to_be_deleted:
      logging.debug('Deleting lost HRU %s', hru)
      band.delete_hru(hru)
    logging.debug('Number of remaining HRUs in band %s is %s', band_id, band.num_hrus)

    # If the band has no current area fraction within the cell, set the
    # median elevation to the floor of the range
    if band.area_frac == 0:
      band.median_elev = band.lower_bound

    logging.debug('Finished state update for band %s.', band_id)

  else:
    logging.debug('No changes in band or glacier area fractions were found for '
      'band %s. No state changes applied.', band_id)
  logging.debug('Number of remaining HRUs in band %s is %s',
                band_id, band.num_hrus)


def update_hru_state(source_hru, dest_hru, case, **kwargs):
  """ Updates the set of state variables for a given HRU based on which of the
    5 cases from the State Update Spec 3.0 is true.
  """
  def carry_over(source, dest):
    dest = source

  if case == '1':
    pass
  elif case == '2':
    pass
  elif case == '3':
    for var in source_hru.hru_state.variables:
      if var in spec_2_vars:
        if type(source_hru.hru_state.variables[var]) == np.ndarray:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == np.ndarray:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / kwargs['new_hru_area_frac'])
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / kwargs['new_hru_area_frac'])
        else:
          dest_hru.hru_state.variables[var] = source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / kwargs['new_hru_area_frac'])
      elif var in spec_3_vars: # SNOW_DENSITY
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0: # avoid division by zero
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars:
        dest_hru.hru_state.variables[var] = source_hru.hru_state.variables[var] \
        * (source_hru.area_frac / kwargs['new_hru_area_frac'])
      elif var in spec_5_vars:
        # TODO: if passed the veg_type, we could set this to zero for non-glacier dest_hrus
        # rather than carrying over a nan
        carry_over(source_hru.hru_state.variables[var], \
          dest_hru.hru_state.variables[var])
      elif var in spec_6_vars: # SNOW_COLD_CONTENT
        dest_hru.hru_state.variables[var] \
        = (dest_hru.hru_state.variables['SNOW_SURF_TEMP'] \
          * (min(MAX_SURFACE_SWE, dest_hru.hru_state.variables['SNOW_SWQ'])) \
          * CH_ICE)
      elif var in spec_7_vars:
        carry_over(source_hru.hru_state.variables[var], \
          dest_hru.hru_state.variables[var])
      elif var in spec_9_vars:
        carry_over(source_hru.hru_state.variables[var], \
          dest_hru.hru_state.variables[var])

  elif case == '4a': # transferring state to open ground
    for var in source_hru.hru_state.variables:
      if var in spec_2_vars:
        if type(source_hru.hru_state.variables[var]) == np.ndarray:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == np.ndarray:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / kwargs['new_open_ground_area_frac'])
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / kwargs['new_open_ground_area_frac'])
        else:
          dest_hru.hru_state.variables[var] \
          = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / kwargs['new_open_ground_area_frac'])
        if var == 'SNOW_CANOPY': # spec-10 sanity check
          if dest_hru.hru_state.variables[var] > 0:
            dest_hru.hru_state.variables['SNOW_SWQ'] += \
              dest_hru.hru_state.variables['SNOW_CANOPY']
            dest_hru.hru_state.variables['SNOW_CANOPY'] = 0
      elif var in spec_3_vars: # SNOW_DENSITY
        # avoid division by zero
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0:
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars: # GLAC_WATER_STORAGE
        dest_hru.hru_state.variables[var] \
        = dest_hru.hru_state.variables[var] \
        + source_hru.hru_state.variables[var] \
        * (source_hru.area_frac / kwargs['new_open_ground_area_frac'])
        if dest_hru.hru_state.variables[var] > 0: # spec-10 sanity check
        #FIXME: hard-coding the dist dim to [0] for now
          dest_hru.hru_state.variables['LAYER_MOIST'][0][Cell.Nlayers-1] += \
            dest_hru.hru_state.variables['GLAC_WATER_STORAGE']
          dest_hru.hru_state.variables['GLAC_WATER_STORAGE'] = 0
      elif var in spec_5_vars: # GLAC_CUM_MASS_BALANCE, only applies to glacier HRUs.
        # Glacier HRU gone, but shadow remains, so we have to set state to zero
        # rather than relying on deletion of the HRU to effectively do this
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0
      elif var in spec_6_vars: # SNOW_COLD_CONTENT
        dest_hru.hru_state.variables[var] \
        = (dest_hru.hru_state.variables['SNOW_SURF_TEMP'] \
          * (min(MAX_SURFACE_SWE, dest_hru.hru_state.variables['SNOW_SWQ'])) \
          * CH_ICE)
      elif var in spec_7_vars:
        dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables[var] \
          * dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.hru_state.variables[var] * source_hru.area_frac \
          * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0))) \
          / (dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.area_frac * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0)))
        if var == 'SNOW_LAST_SNOW' or var == 'SNOW_MELTING':
          dest_hru.hru_state.variables[var] = ceil(dest_hru.hru_state.variables[var])
        source_hru.hru_state.variables[var] = 0
      elif var in spec_9_vars:
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0

  elif case == '4b': # transferring state to glacier
    for var in source_hru.hru_state.variables:
      if var in spec_2_vars:
        if type(source_hru.hru_state.variables[var]) == np.ndarray:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == np.ndarray:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / kwargs['new_glacier_area_frac'])
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx]\
              * (source_hru.area_frac / kwargs['new_glacier_area_frac'])
        else:
          dest_hru.hru_state.variables[var] \
          = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / kwargs['new_glacier_area_frac'])
        if var == 'SNOW_CANOPY': # spec-10 sanity check
          if dest_hru.hru_state.variables[var] > 0:
            dest_hru.hru_state.variables['SNOW_SWQ'] += \
              dest_hru.hru_state.variables['SNOW_CANOPY']
            dest_hru.hru_state.variables['SNOW_CANOPY'] = 0
      elif var in spec_3_vars: # SNOW_DENSITY
        # avoid division by zero
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0:
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars: # GLAC_WATER_STORAGE
        pass # this value is carried over from the previous time step
      elif var in spec_5_vars: # GLAC_CUM_MASS_BALANCE (glacier HRUs only)
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0
      elif var in spec_6_vars: # SNOW_COLD_CONTENT
        dest_hru.hru_state.variables[var] \
        = (dest_hru.hru_state.variables['SNOW_SURF_TEMP']
          * (min(MAX_SURFACE_SWE, dest_hru.hru_state.variables['SNOW_SWQ']))
          * CH_ICE)
      elif var in spec_7_vars:
        dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables[var]
          * dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.hru_state.variables[var] * source_hru.area_frac \
          * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0))) \
          / (dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.area_frac * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0)))
        if var == 'SNOW_LAST_SNOW' or var == 'SNOW_MELTING':
          dest_hru.hru_state.variables[var] = ceil(dest_hru.hru_state.variables[var])
        source_hru.hru_state.variables[var] = 0
      elif var in spec_9_vars:
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0

  elif case == '5a' or case == '5d': # transferring state to glacier
    for var in source_hru.hru_state.variables:
      if var in spec_2_vars:
        if type(source_hru.hru_state.variables[var]) == np.ndarray:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == np.ndarray:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / kwargs['new_glacier_area_frac'])
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / kwargs['new_glacier_area_frac'])
        else:
          dest_hru.hru_state.variables[var] = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / kwargs['new_glacier_area_frac'])
        if var == 'SNOW_CANOPY': # spec-10 sanity check
          if dest_hru.hru_state.variables[var] > 0:
            dest_hru.hru_state.variables['SNOW_SWQ'] += \
              dest_hru.hru_state.variables['SNOW_CANOPY']
            dest_hru.hru_state.variables['SNOW_CANOPY'] = 0
      elif var in spec_3_vars: # SNOW_DENSITY
         # avoid division by zero
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0:
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars:
        dest_hru.hru_state.variables[var] \
        = dest_hru.hru_state.variables[var] \
        + source_hru.hru_state.variables[var] \
        * (source_hru.area_frac / kwargs['new_glacier_area_frac'])
      elif var in spec_5_vars: # GLAC_CUM_MASS_BALANCE (glacier HRUs only)
        # Glacier HRU gone, but shadow remains, so we have to set state to zero
        # rather than relying on deletion of the HRU to effectively do this
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0
      elif var in spec_6_vars: # SNOW_COLD_CONTENT
        dest_hru.hru_state.variables[var] \
        = (dest_hru.hru_state.variables['SNOW_SURF_TEMP'] \
          * (min(MAX_SURFACE_SWE, dest_hru.hru_state.variables['SNOW_SWQ'])) \
          * CH_ICE)
      elif var in spec_7_vars:
        dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables[var] \
          * dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.hru_state.variables[var] * source_hru.area_frac \
          * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0))) \
          / (dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.area_frac * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0)))
        if var == 'SNOW_LAST_SNOW' or var == 'SNOW_MELTING':
          dest_hru.hru_state.variables[var] = ceil(dest_hru.hru_state.variables[var])
        source_hru.hru_state.variables[var] = 0
      elif var in spec_9_vars:
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0

  elif case == '5b': # transferring state to open ground
    for var in source_hru.hru_state.variables:
      if var in spec_2_vars:
        if type(source_hru.hru_state.variables[var]) == np.ndarray:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == np.ndarray:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / kwargs['new_open_ground_area_frac'])
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / kwargs['new_open_ground_area_frac'])
        else:
          dest_hru.hru_state.variables[var] \
          = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / kwargs['new_open_ground_area_frac'])
        if var == 'SNOW_CANOPY': # spec-10 sanity check
          if dest_hru.hru_state.variables[var] > 0:
            dest_hru.hru_state.variables['SNOW_SWQ'] += \
              dest_hru.hru_state.variables['SNOW_CANOPY']
            dest_hru.hru_state.variables['SNOW_CANOPY'] = 0
      elif var in spec_3_vars: # SNOW_DENSITY
         # avoid division by zero
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0:
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars:
        dest_hru.hru_state.variables[var] \
        = dest_hru.hru_state.variables[var] \
        + source_hru.hru_state.variables[var] \
        * (source_hru.area_frac / kwargs['new_open_ground_area_frac'])
        if var == 'GLAC_WATER_STORAGE': # spec-10 sanity check
          if dest_hru.hru_state.variables[var] > 0:
            #FIXME: hard-coding the dist dim to [0] for now
            dest_hru.hru_state.variables['LAYER_MOIST'][0][Cell.Nlayers-1] += \
              dest_hru.hru_state.variables['GLAC_WATER_STORAGE']
            dest_hru.hru_state.variables['GLAC_WATER_STORAGE'] = 0
      elif var in spec_5_vars: # GLAC_CUM_MASS_BALANCE (glacier HRUs only)
        # Glacier HRU gone, but shadow remains, so we have to set state to zero
        # rather than relying on deletion of the HRU to effectively do this
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0
      elif var in spec_6_vars: # SNOW_COLD_CONTENT
        dest_hru.hru_state.variables[var] \
        = (dest_hru.hru_state.variables['SNOW_SURF_TEMP'] \
          * (min(MAX_SURFACE_SWE, dest_hru.hru_state.variables['SNOW_SWQ'])) \
          * CH_ICE)
      elif var in spec_7_vars:
        dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables[var] \
          * dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.hru_state.variables[var] * source_hru.area_frac \
          * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0))) \
          / (dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.area_frac * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0)))
        if var == 'SNOW_LAST_SNOW' or var == 'SNOW_MELTING':
          dest_hru.hru_state.variables[var] = ceil(dest_hru.hru_state.variables[var])
        source_hru.hru_state.variables[var] = 0
      elif var in spec_9_vars:
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0

  elif case == '5c': # transferring state to vegetated HRU
    for var in source_hru.hru_state.variables:
      if var in spec_2_vars:
        if type(source_hru.hru_state.variables[var]) == np.ndarray:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == np.ndarray:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / kwargs['new_hru_area_frac'])
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / kwargs['new_hru_area_frac'])
        else:
          dest_hru.hru_state.variables[var] \
          = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / kwargs['new_hru_area_frac'])
        if var == 'SNOW_CANOPY': # spec-10 sanity check
          if dest_hru.hru_state.variables[var] > 0:
            dest_hru.hru_state.variables['SNOW_SWQ'] += \
              dest_hru.hru_state.variables['SNOW_CANOPY']
            dest_hru.hru_state.variables['SNOW_CANOPY'] = 0
      elif var in spec_3_vars: # SNOW_DENSITY
         # avoid division by zero
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0:
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars:
        dest_hru.hru_state.variables[var] \
        = dest_hru.hru_state.variables[var] \
        + source_hru.hru_state.variables[var] \
        * (source_hru.area_frac / kwargs['new_hru_area_frac'])
        if var == 'GLAC_WATER_STORAGE': # spec-10 sanity check
          if dest_hru.hru_state.variables[var] > 0:
            #FIXME: hard-coding the dist dim to [0] for now
            dest_hru.hru_state.variables['LAYER_MOIST'][0][Cell.Nlayers-1] += \
              dest_hru.hru_state.variables['GLAC_WATER_STORAGE']
            dest_hru.hru_state.variables['GLAC_WATER_STORAGE'] = 0
      elif var in spec_5_vars: # GLAC_CUM_MASS_BALANCE (glacier HRUs only)
        # Glacier HRU gone, but shadow remains, so we have to set state to zero
        # rather than relying on deletion of the HRU to effectively do this
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0
      elif var in spec_6_vars: # SNOW_COLD_CONTENT
        dest_hru.hru_state.variables[var] \
        = (dest_hru.hru_state.variables['SNOW_SURF_TEMP'] \
          * (min(MAX_SURFACE_SWE, dest_hru.hru_state.variables['SNOW_SWQ'])) \
          * CH_ICE)
      elif var in spec_7_vars:
        dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables[var] \
          * dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.hru_state.variables[var] * source_hru.area_frac \
          * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0))) \
          / (dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.area_frac * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0)))
        if var == 'SNOW_LAST_SNOW' or var == 'SNOW_MELTING':
          dest_hru.hru_state.variables[var] = ceil(dest_hru.hru_state.variables[var])
        source_hru.hru_state.variables[var] = 0
      elif var in spec_9_vars:
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0
