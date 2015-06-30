''' cells.py

    This module captures the conceptual components of the VIC grid cell,
    broken into classes by elevation band and Hydrological Response Unit (HRU).
    It also provides functions that operate on these class objects, nested
    within the OrderedDict of VIC cells. 

'''

from collections import OrderedDict
from copy import deepcopy
import numpy as np

class Band(object):
    """ Class capturing VIC cell parameters at the elevation band level
    """
    # glacier_id and open_ground_id, glacier_root_zone_parms, open_ground_root_zone_parms,
    # and band_size are defined on a *per-run* basis.
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
        """ The Band area fraction, equal to the sum of HRU area fractions within this band, 
            which should be equal to the total area fraction of the band as given in the 
            Snow Band Parameters file
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
        return sum([hru.area_frac for veg_type, hru in self.hrus.items() if veg_type == self.open_ground_id])

    def create_hru(self, veg_type, area_frac):
        """ Creates a new HRU of veg_type
        """
        # Append new_hru to existing list of HRUs for this band
        if veg_type == self.glacier_id:
            self.hrus[veg_type] = HydroResponseUnit(area_frac, self.glacier_root_zone_parms)
        elif veg_type == self.open_ground_id:
            self.hrus[veg_type] = HydroResponseUnit(area_frac, self.open_ground_root_zone_parms)

    def delete_hru(self, veg_type):
        """ Deletes an HRU of veg_type within the Band
        """
        del self.hrus[veg_type]

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__, self.median_elev,
                                   repr(self.hrus))
    def __str__(self):
        return '{}(@{} meters with {} HRUs)'.format(self.__class__.__name__,
                                                    self.median_elev,
                                                    len(self.hrus))

def apply_custom_root_zone_parms(hru_cell_dict, glacier_root_zone_parms, open_ground_root_zone_parms):
    """ Utility function to apply user-supplied custom root zone parameters to glacier 
        and/or open ground HRUs at initialization time """
    for cell in hru_cell_dict.keys():
        for key in hru_cell_dict[cell]:
            if glacier_root_zone_parms and (key[1] == global_parms.glacier_id):
                hru_cell_dict[cell][key].root_zone_parms = glacier_root_zone_parms
            if open_ground_root_zone_parms and (key[1] == global_parms.open_ground_id):
                hru_cell_dict[cell][key].root_zone_parms = open_ground_root_zone_parms

def merge_cell_input(hru_cell_dict, elevation_cell_dict):
    """ Utility function to merge the dict of HRUs loaded via vegparams.load_veg_parms()
        with the PaddedDeque of Bands loaded at start-up via snbparams.load_snb_parms()
        into one unified structure capturing all VIC cells' properties
    """ 
    missing_keys = hru_cell_dict.keys() ^ elevation_cell_dict.keys()
    if missing_keys:
        raise Exception("One or more cell IDs were found in one input file,"
                "but not the other. IDs: {}".format(missing_keys))

    # initialize new cell container
    cells = deepcopy(elevation_cell_dict)
    # FIXME: this is a little awkward
    for cell_id, hru_dict in hru_cell_dict.items():
        band_ids = { band_id for band_id, _ in hru_dict.keys() }
        band_dict = { band_id: {} for band_id in band_ids }
        for (band_id, veg_type), hru in hru_dict.items():
            band_dict[band_id][veg_type] = hru
        for band_id, hru_dict in band_dict.items():
            cells[cell_id][band_id].hrus = hru_dict
    return cells

class HydroResponseUnit(object):
    """ Class capturing vegetation parameters at the single vegetation
        tile (HRU) level (of which there can be many per band)
    """
    def __init__(self, area_frac, root_zone_parms):
        self.area_frac = area_frac
        self.root_zone_parms = root_zone_parms
    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__,
                                   self.area_frac, self.root_zone_parms)
    def __str__(self):
        return 'HRU({:.2f}%, {})'.format(self.area_frac * 100, self.root_zone_parms)

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and self.__dict__ == other.__dict__)

def update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,
                      surf_dem, num_rows_dem, num_cols_dem, glacier_mask):
    """ Applies the updated RGM DEM and glacier mask and calculates and updates all HRU area fractions 
        for all elevation bands within the VIC cells
    """
    # Create temporary band area and glacier area bins
    band_areas = {} # count of pixels; a proxy for area, within each band
    glacier_areas = {} # count of pixels landing within the glacier mask; a proxy for glacier area

    # Scan through the RGM output DEM and drop pixels into bins according to cell and elevation
    for cell_id, cell in cells.items():
        band_areas[cell_id] = [0] * num_snow_bands
        glacier_areas[cell_id] = [0] * num_snow_bands

        # Select the portion of the DEM that pertains to this cell from which we will do binning
        masked_dem = np.ma.masked_array(surf_dem)
        masked_dem[np.where(cellid_map != float(cell_id))] = np.ma.masked
        flat_dem = masked_dem[~masked_dem.mask]

        # Check if any pixels fall outside of valid range of bands
        if len(np.where(flat_dem < cell[0].lower_bound)[0]) > 0:
            raise Exception(
                'One or more RGM output DEM pixels lies below the lowest \
                defined elevation band (< {}m) as defined by the Snow Band Parameter File \
                for cell {}. You may need to add or shift the zero padding to accommodate this.'
                ).format(cell[0].lower_bound, cell_id)
        if len(np.where(flat_dem > cell[-1].upper_bound)[0]) > 0:
            raise Exception(
                'One or more RGM output DEM pixels lies above the highest \
                defined elevation band (> {}m) as defined by the Snow Band Parameter File \
                for cell {}. You may need to add or shift the zero padding to accommodate this.'
                ).format(cell[-1].upper_bound, cell_id)

        # Identify the bounds of the band area and glacier area bins for this cell
        # and update the band median elevations 
        band_bin_bounds = []
        for band in cell:
            band_bin_bounds.append(band.lower_bound)
            band_pixels = flat_dem[np.where((flat_dem >= band.lower_bound) & (flat_dem < band.upper_bound))]
            if not band_pixels.size: # if there are no pixels in this band
                band.median_elev = band.lower_bound
            else:
                band.median_elev = np.median(band_pixels)
            print ('band: {}, band.median_elev: {}, band.area_frac: {}, band.area_frac_glacier: {}'.format(band, band.median_elev, band.area_frac, band.area_frac_glacier))
        print('band_bin_bounds: {}'.format(band_bin_bounds))

        # Select portion of glacier_mask pertaining to this cell
        cell_glacier_mask = np.ma.masked_array(glacier_mask)
        cell_glacier_mask[np.where(cellid_map != float(cell_id))] = np.ma.masked
        masked_dem.mask = False
        masked_dem[np.where(cell_glacier_mask !=1)] = np.ma.masked
        flat_glacier_dem = masked_dem[~masked_dem.mask]

        # Do binning into band_areas[cell][band] and glacier_areas[cell][band] 
        # using data in flat_dem and flat_glacier_dem
        inds = np.digitize(flat_dem, band_bin_bounds)
        print('band bin inds-1: {}'.format(inds-1))
        for idx, count in enumerate(np.bincount(inds-1)):
            band_areas[cell_id][idx] = count

        inds = np.digitize(flat_glacier_dem, band_bin_bounds)
        print('glacier bin inds-1: {}'.format(inds-1))
        for idx, count in enumerate(np.bincount(inds-1)):
            glacier_areas[cell_id][idx] = count
        print ('band_areas: {}'.format(band_areas))
        print('glacier_areas: {}'.format(glacier_areas))
        
        # Update all Band and HRU area fractions in this cell
        for band_id, band in enumerate(cell):
            # Update total area fraction for this Band
            new_band_area_frac = band_areas[cell_id][band_id] / cell_areas[cell_id] 
            new_glacier_area_frac = glacier_areas[cell_id][band_id] / cell_areas[cell_id]
            # If the glacier HRU area fraction has changed for this band then we need to update all area fractions
            if new_glacier_area_frac != cell[band_id].area_frac_glacier:
                # Calculate new non-glacier area fraction for this band
                new_non_glacier_area_frac = new_band_area_frac - new_glacier_area_frac
                # Calculate new_residual area fraction
                new_residual_area_frac = new_non_glacier_area_frac - cell[band_id].area_frac_non_glacier
    # NOTE: what should be done here for Bands just created above?  cells[cell][band_idx].area_frac_non_glacier = 0 & cells[cell][band_idx].area_frac_open_ground = 1 ?
                # Calculate new open ground area fraction
                new_open_ground_area_frac = np.max([0, (cell[band_id].area_frac_open_ground + new_residual_area_frac)])
                
                # Use old proportions of vegetated areas for scaling their area fractions in this iteration
                veg_scaling_divisor = cell[band_id].area_frac_non_glacier - cell[band_id].area_frac_open_ground
                # Calculate the change in sum of vegetated area fractions
                delta_area_vegetated = np.min([0, (cell[band_id].area_frac_open_ground + new_residual_area_frac)])

                # Update area fractions for all HRUs in this Band  
                glacier_found = False
                open_ground_found = False
                hrus_to_be_deleted = []
                for veg_type, hru in band.hrus.items():
                    print('cell_id: {}, band_id: {}, veg_type: {}'.format(cell_id, band_id, veg_type))
                    if veg_type == Band.glacier_id:
                        print('glacier found')
                        glacier_found = True
                        cell[band_id].hrus[veg_type].area_frac = new_glacier_area_frac
                    elif veg_type == Band.open_ground_id:
                        open_ground_found = True
                        # If open ground area fraction was reduced to 0, we delete the HRU
                        if new_open_ground_area_frac == 0:
                            hrus_to_be_deleted.append(Band.open_ground_id)
                            print('OPEN GROUND HRU ({}) MARKED FOR DELETION'.format(veg_type))
                            # band.delete_hru(Band.open_ground_id)
                        else:
                            cell[band_id].hrus[veg_type].area_frac = new_open_ground_area_frac
                    else:
                        print('veg section... cell_id: {}, band_id: {}, veg_type: {}'.format(cell_id, band_id, veg_type))
                        # Calculate change in HRU area fraction & update
                        delta_area_hru = delta_area_vegetated * (cell[band_id].hrus[veg_type].area_frac / veg_scaling_divisor)
                        new_hru_area_frac = cell[band_id].hrus[veg_type].area_frac + delta_area_hru
                        if new_hru_area_frac == 0: # HRU has disappeared, never to return (only open ground can come back in its place)
                            hrus_to_be_deleted.append(veg_type)
                            print('VEGETATED HRU ({}) MARKED FOR DELETION'.format(veg_type))
                            # band.delete_hru(veg_type)
                        else:
                            cell[band_id].hrus[veg_type].area_frac = new_hru_area_frac
                
                # Remove HRUs marked for deletion
                for hru in hrus_to_be_deleted:
                    band.delete_hru(hru)
                # Create glacier HRU if it didn't already exist, if we have a non-zero new_glacier_area_frac. 
                # If glacier area fraction was reduced to 0, we leave the "shadow glacier" HRU in place for VIC
                if not glacier_found and (new_glacier_area_frac > 0):
                    band.create_hru(Band.glacier_id, new_glacier_area_frac)
                # If open ground was exposed in this band we need to add an open ground HRU
                if not open_ground_found and (new_open_ground_area_frac > 0):
                    band.create_hru(Band.open_ground_id, new_open_ground_area_frac)
                # If there are no HRUs left in this band, it's now a dummy band (set its median elevation to the floor)
                if band.num_hrus == 0:
                    band.median_elev = band.lower_bound

                # Sanity check that glacier + non-glacier area fractions add up to Band's total area fraction, within tolerance
                sum_test = new_glacier_area_frac + new_non_glacier_area_frac
                if sum_test != new_band_area_frac:
                    raise Exception(
                            'cell {}, band {}: glacier area fraction {} + '
                            'non-glacier area fraction {} = {} is not equal to '
                            'the band area fraction of {}'
                            .format(cell_id, band_id, new_glacier_area_frac,
                                    new_non_glacier_area_frac, sum_test,
                                    new_band_area_frac)
                    )
