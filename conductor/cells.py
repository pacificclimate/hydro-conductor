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

# These are set in the VIC header snow.h, but should probably be passed in via the state file
MAX_SURFACE_SWE = 0.125
CH_ICE = 2100E+03

class Cell(object):
  """Class capturing VIC cells
  """
  def __init__(self, bands):
    self.cell_state = CellState()
    self.bands = bands

  def __eq__(self, other):
    return (self.__class__ == other.__class__ and self.__dict__ == other.__dict__)

class CellState(object):
  """Class capturing the set of VIC cell state and metadata variables that can
    change in a yearly VIC run.
  """
  def __init__(self):
    self.variables = {
      'SOIL_DZ_NODE': [],
      'SOIL_ZSUM_NODE': [],
      'GRID_CELL':-1,
      'NUM_BANDS': 0,
      'VEG_TYPE_NUM': 0,
      'GLAC_MASS_BALANCE_EQN_TERMS': []
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

  def create_hru(self, veg_type, area_frac):
    """Creates a new HRU of veg_type
    """
    # Append new hru to existing dict of HRUs for this band
    if veg_type == self.glacier_id:
      self.hrus[veg_type] = HydroResponseUnit(area_frac,\
        self.glacier_root_zone_parms)
    elif veg_type == self.open_ground_id:
      self.hrus[veg_type] = HydroResponseUnit(area_frac,\
        self.open_ground_root_zone_parms)

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
  def __init__(self, area_frac, root_zone_parms):
    self.area_frac = area_frac
    self.root_zone_parms = root_zone_parms
    self.hru_state = HruState()
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
  def __init__(self):
    # variables is an OrderedDict because there is temporal dependence in the
    # state update among some of them when update_hru_state() is called
    self.variables = OrderedDict([
      # HRU state variables with dimensions (lat, lon, hru)
      ('HRU_BAND_INDEX', -1),
      ('HRU_VEG_INDEX', -1),
      # These two have dimensions (lat, lon, hru, dist, Nlayers)
      ('LAYER_ICE_CONTENT', []),
      ('LAYER_MOIST', []),
      # HRU_VEG_VAR_WDEW has dimensions (lat, lon, hru, dist)
      ('HRU_VEG_VAR_WDEW' , []),
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
      ('ENERGY_T', []),
      ('ENERGY_T_FBCOUNT', []),
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
spec_1_vars = ['NUM_BANDS', 'SOIL_DZ_NODE', 'SOIL_ZSUM_NODE',\
  'VEG_TYPE_NUM', 'HRU_BAND_INDEX', 'HRU_VEG_INDEX']

spec_2_vars = ['LAYER_ICE_CONTENT', 'LAYER_MOIST',\
  'HRU_VEG_VAR_WDEW', 'SNOW_CANOPY', 'SNOW_DEPTH',\
  'SNOW_PACK_WATER', 'SNOW_SURF_WATER', 'SNOW_SWQ',\
  'SNOW_PACK_TEMP', 'SNOW_SURF_TEMP']

spec_3_vars = ['SNOW_DENSITY']

spec_4_vars = ['GLAC_WATER_STORAGE']

spec_5_vars = ['GLAC_CUM_MASS_BALANCE']

spec_6_vars = ['SNOW_COLD_CONTENT']

spec_7_vars = ['SNOW_ALBEDO', 'SNOW_LAST_SNOW', 'SNOW_MELTING']

# spec_8_vars = ['SNOW_LAST_SNOW', 'SNOW_MELTING']

spec_9_vars = ['ENERGY_T', 'ENERGY_TFOLIAGE', 'GLAC_SURF_TEMP',\
  'ENERGY_TCANOPY_FBCOUNT', 'ENERGY_TFOLIAGE_FBCOUNT','ENERGY_TSURF_FBCOUNT',\
  'GLAC_SURF_TEMP_FBCOUNT', 'SNOW_SURF_TEMP_FBCOUNT',\
  # miscellaneous vars:
  'GLAC_SURF_TEMP_FBFLAG', 'GLAC_VAPOR_FLUX',\
  'SNOW_CANOPY_ALBEDO', 'SNOW_SURFACE_FLUX', 'SNOW_SURF_TEMP_FBFLAG',\
  'SNOW_TMP_INT_STORAGE', 'SNOW_VAPOR_FLUX']
# NOTE: MAY STILL NEED TO ADD DEFERRED VARS TO SPEC_9_VARS

spec_10_vars = ['SNOW_CANOPY', 'SNOW_SWQ', 'GLAC_WATER_STORAGE']

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
  return cells

def update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
  surf_dem, num_rows_dem, num_cols_dem, glacier_mask):
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
    band_areas[cell_id] = [0] * num_snow_bands
    glacier_areas[cell_id] = [0] * num_snow_bands

    # Select the portion of the DEM that pertains to this cell from which we
    # will do binning
    masked_dem = np.ma.masked_array(surf_dem)
    masked_dem[np.where(cellid_map != float(cell_id))] = np.ma.masked
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
    cell_glacier_mask[np.where(cellid_map != float(cell_id))] = np.ma.masked
    masked_dem.mask = False
    masked_dem[np.where(cell_glacier_mask !=1)] = np.ma.masked
    # Create a regular 'flat' np.array of the DEM subset that is the glacier mask
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
    previous_glacier_area_fracs = [0] * len(cell.bands)
    new_band_area_frac = [0] * len(cell.bands)
    new_glacier_area_frac = [0] * len(cell.bands)
    new_non_glacier_area_frac = [0] * len(cell.bands)
    new_open_ground_area_frac = [0] * len(cell.bands)
    veg_scaling_divisor = [0] * len(cell.bands)
    delta_area_vegetated = [0] * len(cell.bands)
    new_hru_area_frac = {}

    ### Update all Band area fractions for this cell
    # This must be done for all bands before the HRU area
    # fractions and state update steps
    for band_id, band in reverse_enumerate(cell.bands):
      # Update total area fraction for this Band
      new_band_area_frac[band_id] = band_areas[cell_id][band_id]/cell_areas[cell_id]
      new_glacier_area_frac[band_id] = glacier_areas[cell_id][band_id]/cell_areas[cell_id]

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
            delta_area_hru = delta_area_vegetated[band_id] * \
              (band.hrus[veg_type].area_frac / veg_scaling_divisor[band_id])
            new_hru_area_frac[str(band_id)] = {}
            new_hru_area_frac[str(band_id)][str(veg_type)] = band.hrus[veg_type].area_frac + delta_area_hru

    ### Update all HRU states for each band, then apply new HRU area
    # fractions calculated above.
    # HRU area fractions from the previous time step are held in the
    # Band / HRU members of a given Cell instance
    # (e.g. band.area_frac_glacier or band.hrus[veg_type].area_frac)
    # and are compared against their newly calculated values for the
    # current time step (e.g. new_glacier_area_frac[band_id]) to determine
    # whether HRUs should be added/deleted, and which state update case
    # should be carried out by update_hru_state().
    # The algorithm requires that we update glacier and open ground HRU
    # states first, followed by all remaining HRU vegetation types.
    for band_id, band in reverse_enumerate(cell.bands):
      hrus_to_be_deleted = [] # HRUs marked for deletion at end of band iteration
      # We need to update area fractions and update state in a band if
      # the band's total area fraction has changed, or the band's glacier
      # HRU area fraction has changed.
      if (new_glacier_area_frac[band_id] != band.area_frac_glacier)\
        or (new_band_area_frac[band_id] != band.area_frac):

        ## Glacier HRU update:
        # NOTE: we never mark glacier HRUs for deletion, we just set their
        # area fractions to zero. This is VIC's 'shadow glacier' requirement.
        if new_glacier_area_frac[band_id] > 0 and band.area_frac_glacier == 0:
          # CASE 1: A new glacier HRU has appeared in this band. Create a new
          # glacier HRU (HRU state defaults are automatically set upon HRU
          # creation, so technically there's no need to call update_hru_state())
          band.create_hru(Band.glacier_id, new_glacier_area_frac[band_id])
          new_area_fracs = {}
          update_hru_state(None,None,'1',**new_area_fracs)
        elif new_glacier_area_frac[band_id] == band.area_frac_glacier:
          # CASE 2: Trivial case. There was no change in glacier HRU area,
          # so we don't call update_hru_state())
          pass
        elif new_glacier_area_frac[band_id] > 0 and band.area_frac_glacier > 0 \
          and new_glacier_area_frac[band_id] != band.area_frac_glacier:
          # CASE 3: Glacier HRU exists in previous and current time steps, but
          # its area fraction has changed.
          new_area_fracs = {
            'new_hru_area_frac': new_glacier_area_frac[band_id]
          }
          update_hru_state( \
            band.hrus[Band.glacier_id], \
            band.hrus[Band.glacier_id], \
            '3', **new_area_fracs)
        elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
          and new_band_area_frac[band_id] > 0 and band.area_frac > 0:
          # CASE 4a: Glacier HRU has disappeared, but the band remains
          # (implies open ground HRU is expanding). Add state to open ground
          # HRU (leave the zero area "shadow glacier" HRU in place for VIC).
          new_area_fracs = {
            'new_open_ground_area_frac': new_open_ground_area_frac[band_id]
          }
          update_hru_state( \
            band.hrus[Band.glacier_id], \
            band.hrus[Band.open_ground_id], \
            '4a', **new_area_fracs)
        elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
          and new_band_area_frac[band_id] == 0 and band.area_frac != 0 \
          and (band_id - 1) >= 0 \
          and new_glacier_area_frac[band_id - 1] > 0:
          # CASE 5a: Glacier HRU and the band have disappeared.
          # An adjacent lower band exists and there's glacier in it
          # (with non-zero area fraction, i.e. not a shadow glacier). Add
          # state to the glacier HRU in the lower band.
          if not Band.glacier_id in cell.bands[band_id - 1].hrus:
            # Due to the top-to-bottom iteration ordering of band processing,
            # we need to create a glacier HRU in the lower band if one doesn't
            # already exist
            cell.bands[band_id - 1].create_hru(Band.glacier_id, new_glacier_area_frac[band_id - 1])
          new_area_fracs = {
            'new_glacier_area_frac': new_glacier_area_frac[band_id - 1]
          }
          update_hru_state( \
            band.hrus[Band.glacier_id], \
            cell.bands[band_id - 1].hrus[Band.glacier_id], \
            '5a', **new_area_fracs)
        elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
          and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
          and (band_id - 1) >= 0 \
          and new_open_ground_area_frac[band_id - 1] > 0:
          # CASE 5b: Glacier HRU and the band have disappeared.
          # An adjacent lower band exists with a non-zero area fraction open
          # ground HRU in it. Add state to that lower band's open ground HRU.
          if not Band.open_ground_id in cell.bands[band_id - 1].hrus:
            # Due to the top-to-bottom iteration ordering of band processing,
            # we need to create an open ground HRU in the lower band if one
            # doesn't already exist
            cell.bands[band_id - 1].create_hru(Band.open_ground_id, new_open_ground_area_frac[band_id - 1])
          new_area_fracs = {
            'new_open_ground_area_frac': new_open_ground_area_frac[band_id - 1]
          }
          update_hru_state( \
            band.hrus[Band.glacier_id], \
            cell.bands[band_id - 1].hrus[Band.open_ground_id], \
            '5b', **new_area_fracs)
        elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0 \
          and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
          and (band_id - 1) >= 0 and new_band_area_frac[band_id - 1] > 0:
          # CASE 5c: Glacier HRU and the band have disappeared due to glacier.
          # If an adjacent lower band exists (with non-zero area fraction), add
          # state to the vegetated HRU with non-zero area fraction and the
          # greatest vegetation type index in that band.
          valid_vegetated_hrus_below = list(filter(lambda x: x[1].area_frac > 0, new_hru_area_frac[str(band_id-1)].items()))
          max_veg_type_below = max(valid_vegetated_hrus_below)[0]
          new_area_fracs = {
            'new_hru_area_frac': new_hru_area_frac[str(band_id - 1)][str(max_veg_type_below)]
          }
          update_hru_state( \
            band.hrus[Band.glacier_id], \
            cell.bands[band_id - 1].hrus[int(max_veg_type_below)], \
            '5c', **new_area_fracs)
        elif new_glacier_area_frac[band_id] == 0 and band.area_frac_glacier > 0\
          and new_band_area_frac[band_id] == 0 and band.area_frac != 0\
          and ((band_id - 1) < 0 or new_band_area_frac[band_id - 1] == 0):
          # CASE 5d: Both the glacier HRU and band have disappeared due
          # to glacier. No adjacent lower band exists (or not one has non-zero
          # area fraction) to transfer state to. In this special case, add
          # state to glacier in the band above.
          new_area_fracs = {
            'new_glacier_area_frac': new_glacier_area_frac[band_id + 1]
          }
          update_hru_state( \
            band.hrus[Band.glacier_id], \
            cell.bands[band_id + 1].hrus[Band.glacier_id], \
            '5d', **new_area_fracs)
        else:
          raise Exception(
            'Error: No state update case identified for cell {}, band {}, HRU {}'
            )

        # Save the glacier area fraction for this band for checking against
        # in the next (lower) band iteration (CASE 5d specifically)
        previous_glacier_area_fracs[band_id] = band.area_frac_glacier

        # Apply update to glacier HRU area fraction
        if Band.glacier_id in band.hrus:
          band.hrus[Band.glacier_id].area_frac = new_glacier_area_frac[band_id]

        ## Open ground HRU update:
        if new_open_ground_area_frac[band_id] > 0 and band.area_frac_open_ground == 0:
          # CASE 1: New open ground was exposed in this band, so
          # create a new open ground HRU (HRU state defaults are
          # automatically set upon HRU creation, so technically
          # there's no need to call update_hru_state())
          band.create_hru(Band.open_ground_id, new_open_ground_area_frac[band_id])
          new_area_fracs = {}
          update_hru_state(None,None,'1',**new_area_fracs)
        elif new_open_ground_area_frac[band_id] == band.area_frac_open_ground:
          # CASE 2: Trivial case. There was no change in open ground HRU area,
          # so we don't call update_hru_state())
          pass
        elif new_open_ground_area_frac[band_id] > 0 and band.area_frac_open_ground > 0 \
          and new_open_ground_area_frac[band_id] != band.area_frac_open_ground:
          # CASE 3: Open ground HRU exists in previous and current time
          # steps, but its area fraction has changed.
          new_area_fracs = {
            'new_hru_area_frac': new_open_ground_area_frac[band_id]
          }
          update_hru_state( \
            band.hrus[Band.open_ground_id], \
            band.hrus[Band.open_ground_id], \
            '3', **new_area_fracs)
        elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
          and new_band_area_frac[band_id] != 0 and band.area_frac > 0:
          # CASE 4b: Open ground has disappeared, but the band remains
          # (implies glacier HRU expanding). Add state to the glacier HRU.
          # Mark the open ground HRU for deletion.
          hrus_to_be_deleted.append(Band.open_ground_id)
          new_area_fracs = {
            'new_glacier_area_frac': new_glacier_area_frac[band_id]
          }
          update_hru_state( \
            band.hrus[Band.open_ground_id], \
            band.hrus[Band.glacier_id], \
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
          hrus_to_be_deleted.append(Band.open_ground_id)
          if not Band.glacier_id in cell.bands[band_id - 1].hrus:
            # Due to the top-to-bottom iteration ordering of band processing,
            # we need to create a glacier HRU in the lower band if one doesn't
            # already exist
            cell.bands[band_id - 1].create_hru(Band.glacier_id, new_glacier_area_frac[band_id - 1])
          new_area_fracs = {
            'new_glacier_area_frac': new_glacier_area_frac[band_id - 1]
          }
          update_hru_state( \
            band.hrus[Band.open_ground_id], \
            cell.bands[band_id - 1].hrus[Band.glacier_id], \
            '5a', **new_area_fracs)
        elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
          and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
          and (band_id - 1) >= 0 \
          and new_open_ground_area_frac[band_id - 1] > 0:
          # CASE 5b: Open ground HRU and the band have disappeared.
          # An adjacent lower band exists with a non-zero area fraction open
          # ground HRU in it. Add state to that lower band's open ground HRU.
          # Mark the open ground HRU for deletion.
          if not Band.open_ground_id in cell.bands[band_id - 1].hrus:
            # Due to the top-to-bottom iteration ordering of band processing,
            # we need to create an open ground HRU in the lower band if one
            # doesn't already exist
            cell.bands[band_id - 1].create_hru(Band.open_ground_id, new_open_ground_area_frac[band_id - 1])
          hrus_to_be_deleted.append(Band.open_ground_id)
          new_area_fracs = {
            'new_open_ground_area_frac': new_open_ground_area_frac[band_id - 1]
          }
          update_hru_state( \
            band.hrus[Band.open_ground_id], \
            cell.bands[band_id - 1].hrus[Band.open_ground_id], \
            '5b', **new_area_fracs)
        elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
          and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
          and (band_id - 1) >= 0 and new_band_area_frac[band_id - 1] > 0:
          # CASE 5c: Open ground HRU and the band have disappeared due to
          # glacier. Mark the open ground HRU for deletion. If an adjacent
          # lower band exists (with non-zero area_frac), add state to the
          # vegetated HRU with non-zero area_frac and the greatest vegetation
          # type index in that band.  
          hrus_to_be_deleted.append(Band.open_ground_id)
          valid_vegetated_hrus_below = list(filter(lambda x: x[1].area_frac > 0, new_hru_area_frac[str(band_id-1)].items()))
          max_veg_type_below = max(valid_vegetated_hrus_below)[0]
          new_area_fracs = {
            'new_hru_area_frac': new_hru_area_frac[str(band_id - 1)][str(max_veg_type_below)]
          }
          update_hru_state( \
            band.hrus[Band.open_ground_id], \
            cell.bands[band_id - 1].hrus[int(max_veg_type_below)], \
            '5c', **new_area_fracs)
        elif new_open_ground_area_frac[band_id] == 0 and band.area_frac_open_ground > 0 \
          and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
          and ((band_id - 1) < 0 or new_band_area_frac[band_id - 1] == 0):
          # CASE 5d: Both the open ground HRU and band have disappeared due
          # to glacier. Mark the HRU for deletion. No adjacent lower band
          # exists (or not one has non-zero area fraction) to transfer
          # state to. In this special case, add state to glacier in the
          # band above.
          hrus_to_be_deleted.append(Band.open_ground_id)
          new_area_fracs = {
            'new_glacier_area_frac': new_glacier_area_frac[band_id + 1]
          }
          update_hru_state( \
            band.hrus[Band.open_ground_id], \
            cell.bands[band_id + 1].hrus[Band.glacier_id], \
            '5d', **new_area_fracs)
        else:
          raise Exception(
            'Error: No state update case identified for cell {}, band {}, HRU {}'
            )
        # Apply update to open ground HRU area fraction. HRU will get deleted later
        # if this is 0 (it will have already been added to hrus_to_be_deleted)
        if Band.open_ground_id in band.hrus:
          band.hrus[Band.open_ground_id].area_frac = new_open_ground_area_frac[band_id]

        ### Update area fractions and states for all remaining non-glacier and
        # non-open ground HRUs in this band 
        for veg_type, hru in band.hrus.items():
          if veg_type is not Band.glacier_id and veg_type is not Band.open_ground_id:
            if new_hru_area_frac[str(band_id)][str(veg_type)] < 0:
              raise Exception(
              'Error: cell {}, band {}: HRU {} has a negative area fraction ({})'
              .format(cell_id, band_id, veg_type, new_hru_area_frac[str(band_id)][str(veg_type)])
              )
            if new_hru_area_frac[str(band_id)][str(veg_type)] == band.hrus[veg_type].area_frac:
              # CASE 2: Trivial case. There was no change in this HRU's area,
              # so we don't call update_hru_state())
              pass
            elif new_hru_area_frac[str(band_id)][str(veg_type)] > 0 and delta_area_hru != 0:
              # CASE 3: Vegetated HRU exists in previous and current time
              # steps, but its area fraction has changed.
              new_area_fracs = {
                'new_hru_area_frac': new_hru_area_frac[str(band_id)][str(veg_type)]
              }
              update_hru_state( \
                hru, \
                hru, \
                '3', **new_area_fracs)
            elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 \
              and new_band_area_frac[band_id] > 0 and band.area_frac > 0:
              # CASE 4b: Vegetated HRU has disappeared, but the band remains
              # (implies glacier HRU expanding). Mark the HRU for deletion
              # (never to return -- only open ground can come back in its place)
              # and add state to the glacier HRU.
              hrus_to_be_deleted.append(veg_type)
              new_area_fracs = {
                'new_glacier_area_frac': new_glacier_area_frac[band_id]
              }
              update_hru_state( \
                hru, \
                band.hrus[Band.glacier_id], \
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
              hrus_to_be_deleted.append(veg_type)
              if not Band.glacier_id in cell.bands[band_id - 1].hrus:
                # Due to the top-to-bottom iteration ordering of band processing,
                # we need to create a glacier HRU in the lower band if one doesn't
                # already exist
                cell.bands[band_id - 1].create_hru(Band.glacier_id, new_glacier_area_frac[band_id - 1])
              new_area_fracs = {
                'new_glacier_area_frac': new_glacier_area_frac[band_id - 1]
              }
              update_hru_state( \
                hru, \
                cell.bands[band_id - 1].hrus[Band.glacier_id], \
                '5a', **new_area_fracs)
            elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 and band.area_frac > 0 \
              and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
              and (band_id - 1) >= 0 \
              and new_open_ground_area_frac[band_id - 1] > 0:
              # CASE 5b: Both the vegetated HRU and the band have disappeared.
              # An adjacent lower band exists with a non-zero area fraction open
              # ground HRU in it. Add state to that lower band's open ground HRU.
              # Mark the vegetated HRU for deletion.
              hrus_to_be_deleted.append(veg_type)
              if not Band.open_ground_id in cell.bands[band_id - 1].hrus:
                # Due to the top-to-bottom iteration ordering of band processing,
                # we need to create an open ground HRU in the lower band if one
                # doesn't already exist
                cell.bands[band_id - 1].create_hru(Band.open_ground_id, new_open_ground_area_frac[band_id - 1])
              new_area_fracs = {
                'new_open_ground_area_frac': new_open_ground_area_frac[band_id - 1]
              }
              update_hru_state( \
                hru, \
                cell.bands[band_id - 1].hrus[Band.open_ground_id], \
                '5b', **new_area_fracs)
            elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 \
              and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
              and (band_id - 1) >= 0 and new_band_area_frac[band_id - 1] > 0:
              # CASE 5c: Both the vegetated HRU and the band have disappeared due to
              # glacier. Mark the vegetated HRU for deletion. An adjacent
              # lower band exists (with non-zero area fraction), where we will
              # add state to the vegetated HRU with non-zero area fraction and
              # the greatest vegetation type index in that band.  
              hrus_to_be_deleted.append(veg_type)
              valid_vegetated_hrus_below = list(filter(lambda x: x[1].area_frac > 0, new_hru_area_frac[str(band_id-1)].items()))
              max_veg_type_below = max(valid_vegetated_hrus_below)[0]
              new_area_fracs = {
                'new_hru_area_frac': new_hru_area_frac[str(band_id - 1)][str(max_veg_type_below)]
              }
              update_hru_state( \
                hru, \
                cell.bands[band_id - 1].hrus[int(max_veg_type_below)], \
                '5c', **new_area_fracs)
            elif new_hru_area_frac[str(band_id)][str(veg_type)] == 0 \
              and new_band_area_frac[band_id] == 0 and band.area_frac > 0 \
              and ((band_id - 1) < 0 or new_band_area_frac[band_id - 1] == 0):
              # CASE 5d: Both the vegetated HRU and band have disappeared due
              # to glacier. Mark the vegetated HRU for deletion. No adjacent
              # lower band exists (or not one has non-zero area fraction) to
              # transfer state to. In this special case, add state to glacier
              # in the band above. 
              hrus_to_be_deleted.append(veg_type)
              new_area_fracs = {
                'new_glacier_area_frac': new_glacier_area_frac[band_id + 1]
              }
              update_hru_state( \
                hru, \
                cell.bands[band_id + 1].hrus[Band.glacier_id], \
                '5d', **new_area_fracs)
            else:
              raise Exception(
                'Error: No state update case identified for cell {}, band {}, HRU {}'
                )
            # Apply update to the HRU's area fraction. HRU will get deleted later
            # if this is 0 (it will have already been added to hrus_to_be_deleted)
            band.hrus[veg_type].area_frac = new_hru_area_frac[str(band_id)][str(veg_type)]

        # Remove HRUs marked for deletion
        for hru in hrus_to_be_deleted:
          band.delete_hru(hru)

        # If the band has no current area fraction within the cell, set the
        # median elevation to the floor of the range
        if band.area_frac == 0:
          band.median_elev = band.lower_bound

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
    print('update_hru_state: case 3')
    for var in source_hru.hru_state.variables:
      if var in spec_1_vars:
        carry_over(source_hru.hru_state.variables[var], \
          dest_hru.hru_state.variables[var])
      elif var in spec_2_vars:
        if type(source_hru.hru_state.variables[var]) == list:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == list:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / new_hru_area_frac)
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / new_hru_area_frac)
        else:
          dest_hru.hru_state.variables[var] = source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / new_hru_area_frac)
      elif var in spec_3_vars: # SNOW_DENSITY
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0: # avoid division by zero
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars:
        dest_hru.hru_state.variables[var] = source_hru.hru_state.variables[var] \
        * (source_hru.area_frac / new_hru_area_frac)
      elif var in spec_5_vars:
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
  elif case == '4a':
    print('update_hru_state: case 4a')
    for var in source_hru.hru_state.variables:
      if var in spec_1_vars:
        pass # need to handle cell metadata state vars differently than HRU state vars
      elif var in spec_2_vars: # transferring state to open ground
        if type(source_hru.hru_state.variables[var]) == list:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == list:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / new_open_ground_area_frac)
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / new_open_ground_area_frac)
        else:
          dest_hru.hru_state.variables[var] \
          = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / new_open_ground_area_frac)
      elif var in spec_3_vars: # SNOW_DENSITY
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0: # avoid division by zero
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars: # transferring state to open ground
        dest_hru.hru_state.variables[var] = dest_hru.hru_state.variables[var] \
        + source_hru.hru_state.variables[var] \
        * (source_hru.area_frac / new_open_ground_area_frac)
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
      elif var in spec_7_vars: # transferring state to open ground
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

  elif case == '4b':
    print('update_hru_state: case 4b')
    for var in source_hru.hru_state.variables:
      if var in spec_1_vars:
        pass # need to handle cell metadata state vars differently than HRU state vars
      elif var in spec_2_vars: # transferring state to glacier
        if type(source_hru.hru_state.variables[var]) == list:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == list:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / new_glacier_area_frac)
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx]\
              * (source_hru.area_frac / new_glacier_area_frac)
        else:
          dest_hru.hru_state.variables[var] \
          = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / new_glacier_area_frac)
      elif var in spec_3_vars: # SNOW_DENSITY
        # avoid division by zero
        if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0:
          dest_hru.hru_state.variables[var] \
          = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) \
          / dest_hru.hru_state.variables['SNOW_DEPTH']
        else:
          dest_hru.hru_state.variables[var] = 0
      elif var in spec_4_vars: # transferring state to glacier
        dest_hru.hru_state.variables[var] = source_hru.hru_state.variables[var]
      elif var in spec_5_vars: # GLAC_CUM_MASS_BALANCE (glacier HRUs only)
        source_hru.hru_state.variables[var] = 0
        dest_hru.hru_state.variables[var] = 0
      elif var in spec_6_vars: # SNOW_COLD_CONTENT
        dest_hru.hru_state.variables[var] \
        = (dest_hru.hru_state.variables['SNOW_SURF_TEMP'] \
          * (min(MAX_SURFACE_SWE, dest_hru.hru_state.variables['SNOW_SWQ'])) \
          * CH_ICE)
      elif var in spec_7_vars: # transferring state to glacier
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
  elif case == '5a':
    print('update_hru_state: case 5a')
    for var in source_hru.hru_state.variables:
      if var in spec_1_vars:
        pass # need to handle cell metadata state vars differently than HRU state vars
      elif var in spec_2_vars: # transferring state to glacier
        if type(source_hru.hru_state.variables[var]) == list:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == list:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / new_glacier_area_frac)
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / new_glacier_area_frac)
        else:
          dest_hru.hru_state.variables[var] = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / new_glacier_area_frac)
      # elif var in spec_3_vars: # SNOW_DENSITY
      #   if dest_hru.hru_state.variables['SNOW_DEPTH'] > 0: # avoid division by zero
      #     dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables['SNOW_SWQ'] * 1000) / dest_hru.hru_state.variables['SNOW_DEPTH']
      #   else:
      #     dest_hru.hru_state.variables[var] = 0
      # elif var in spec_4_vars: # transferring state to glacier
      #   dest_hru.hru_state.variables[var] = dest_hru.hru_state.variables[var] + source_hru.hru_state.variables[var] * (source_hru.area_frac / new_open_ground_area_frac)
      # elif var in spec_5_vars: # GLAC_CUM_MASS_BALANCE, only applies to glacier HRUs
      #   # Glacier HRU gone, but shadow remains, so we have to set state to zero
      #   # rather than relying on deletion of the HRU to effectively do this
      #   source_hru.hru_state.variables[var] = 0
      #   dest_hru.hru_state.variables[var] = 0
      # elif var in spec_6_vars: # SNOW_COLD_CONTENT
      #   dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables['SNOW_SURF_TEMP'] * (min(MAX_SURFACE_SWE, dest_hru.hru_state.variables['SNOW_SWQ'])) * CH_ICE)
      elif var in spec_7_vars: # transferring state to glacier
        dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables[var] \
          * dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.hru_state.variables[var] * source_hru.area_frac \
          * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0))) \
          / (dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.area_frac * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0)))
        if var == 'SNOW_LAST_SNOW' or var == 'SNOW_MELTING':
          dest_hru.hru_state.variables[var] = ceil(dest_hru.hru_state.variables[var])
        source_hru.hru_state.variables[var] = 0

      # elif var in spec_9_vars:
      #   source_hru.hru_state.variables[var] = 0
      #   dest_hru.hru_state.variables[var] = 0
  elif case == '5b':
    print('update_hru_state: case 5b')
    for var in source_hru.hru_state.variables:
      if var in spec_1_vars:
        pass # need to handle cell metadata state vars differently than HRU state vars
      elif var in spec_2_vars: # transferring state to open ground
        if type(source_hru.hru_state.variables[var]) == list:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == list:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / new_open_ground_area_frac)
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / new_open_ground_area_frac)
        else:
          dest_hru.hru_state.variables[var] \
          = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / new_open_ground_area_frac)
        #
        #
      elif var in spec_7_vars: # transferring state to open ground
        dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables[var] \
          * dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.hru_state.variables[var] * source_hru.area_frac \
          * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0))) \
          / (dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.area_frac * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0)))
        if var == 'SNOW_LAST_SNOW' or var == 'SNOW_MELTING':
          dest_hru.hru_state.variables[var] = ceil(dest_hru.hru_state.variables[var])
        source_hru.hru_state.variables[var] = 0
  elif case == '5c':
    print('update_hru_state: case 5c')
    for var in source_hru.hru_state.variables:
      if var in spec_1_vars:
        pass # need to handle cell metadata state vars differently than HRU state vars
      elif var in spec_2_vars: # transferring state to glacier
        if type(source_hru.hru_state.variables[var]) == list:
          for layer_idx, layer in enumerate(source_hru.hru_state.variables[var]):
            if type(source_hru.hru_state.variables[var][layer_idx]) == list:
              for item_idx, item in enumerate(source_hru.hru_state.variables[var][layer_idx]):
                dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                = dest_hru.hru_state.variables[var][layer_idx][item_idx] \
                + source_hru.hru_state.variables[var][layer_idx][item_idx] \
                * (source_hru.area_frac / new_hru_area_frac)
            else:
              dest_hru.hru_state.variables[var][layer_idx] \
              = dest_hru.hru_state.variables[var][layer_idx] \
              + source_hru.hru_state.variables[var][layer_idx] \
              * (source_hru.area_frac / new_hru_area_frac)
        else:
          dest_hru.hru_state.variables[var] \
          = dest_hru.hru_state.variables[var] \
          + source_hru.hru_state.variables[var] \
          * (source_hru.area_frac / new_hru_area_frac)
        #
        #
        #
      elif var in spec_7_vars: # transferring state to vegetated HRU
        dest_hru.hru_state.variables[var] = (dest_hru.hru_state.variables[var] \
          * dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.hru_state.variables[var] * source_hru.area_frac \
          * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0))) \
          / (dest_hru.area_frac * (int(dest_hru.hru_state.variables['SNOW_SWQ'] >= 0)) \
          + source_hru.area_frac * (int(source_hru.hru_state.variables['SNOW_SWQ'] >= 0)))
        if var == 'SNOW_LAST_SNOW' or var == 'SNOW_MELTING':
          dest_hru.hru_state.variables[var] = ceil(dest_hru.hru_state.variables[var])
        source_hru.hru_state.variables[var] = 0
  elif case == '5d':
    print('update_hru_state: case 5d')

