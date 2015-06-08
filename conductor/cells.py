''' cells.py

    This module captures the conceptual components of the VIC grid cell,
    broken into classes by elevation band and Hydrological Response Unit (HRU).
    It also provides functions that operate on these class objects, nested
    within the OrderedDict of VIC cells. 

'''

class Band(object):
    """ Class capturing vegetation parameters at the elevation band level
    """
    # glacier_id and open_ground_id are defined on a *per-run* basis. On one
    # hand, we want to contain these magic numbers within this class, but OTOH
    # it's more overhead than we would want to maintain the value for *every
    # instance* We'll define them statically in the class, and allow
    # applications to set them in the *unlikely* circumstances with which they
    # would be changed E.g.
    # Band.glacier_id = 10
    # my_band_instance = Band()
    # my_band_instance.area_frac_glacier # uses 10 for the comparison value
    glacier_id = 22
    open_ground_id = 19

    def __init__(self, median_elev, hrus={}):
        self.median_elev = median_elev
        self.hrus = hrus

    @property
    def num_hrus(self):
        return len(self.hrus)

    @property
    def area_frac(self): 
        """ The Band area fraction, equal to the sum of HRU area fractions within this band, 
            which should be equal to the total area fraction of the band as given in the 
            Snow Band Parameters file
        """
        return sum([hru.area_frac for hru in self.hrus.values()])

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

    def create_hru(self, veg_type, area_frac, root_zone_parms):
        """ Creates a new HRU of veg_type
        """
        # Append new_hru to existing list of HRUs for this band
        self.hrus[veg_type] = HydroResponseUnit(area_frac, root_zone_parms)

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

class Cell(object):

    def create_band(self, elevation):
        new_band = Band(elevation)
        if elevation < self.bands.peekleft():
            self.bands.appendleft(new_band)
        elif elevation > self.bands.peekright():
            self.bands.append(new_band)
        else:
            raise ValueError("Cannot create a new band of elevation {} since "
                    "bands already exist for the interval {}-{}"
                    .format(elevation, self.bands.peekleft(), self.bands.peekright()))

    def delete_band(self):
        raise NotImplementedError("FIXME: What's the best API for this method?!")

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

# FIXME: update usage of pixel_to_cell_map
def update_area_fracs(cells, cell_areas, num_snow_bands, band_size, band_map, pixel_to_cell_map,
                      surf_dem, num_rows_dem, num_cols_dem, glacier_mask, glacier_id, open_ground_id):
    """ Calculates and updates the area fractions of elevation bands within VIC cells, and
        area fraction of glacier and open ground within VIC cells (broken down by elevation band).
        Calls VegParams::update_tiles() to update all vegetation parameter tiles (or add/delete them as needed)
    """
    all_pixel_elevs = {} # temporarily store pixel elevations in bins by band so the median can be calculated
    band_areas = {} # temporary count of pixels, a proxy for area, within each band
    glacier_areas = {} # temporary count of pixels landing within the glacier mask, a proxy for glacier area
    for cell in cells:
        all_pixel_elevs[cell] = [] * num_snow_bands
        band_areas[cell] = [0] * num_snow_bands
        glacier_areas[cell] = [0] * num_snow_bands
    for row in range(num_rows_dem):
        for col in range(num_cols_dem):
            cell = pixel_to_cell_map[row][col][0] # get the VIC cell this pixel belongs to
            if cell != 'NA':
                # Use the RGM DEM output to update the pixel median elevation in the pixel_to_cell_map
                pixel_elev = float(surf_dem[row][col])
                pixel_to_cell_map[row][col][1] = pixel_elev
                for band_idx, band in enumerate(band_map[cell]):
                    band_found = False
                    if (band < pixel_elev) and (pixel_elev < (band + band_size)):
                        band_found = True
                        # Gather all pixel_elev values to update the median_elev of this band later
                        all_pixel_elevs[cell][band_idx].append(pixel_elev)
                        band_areas[cell][band_idx] += 1
                        if glacier_mask[row][col]:
                            glacier_areas[cell][band_idx] += 1
                        break
                if not band_found: # we have to introduce a new elevation band
                    new_band_idx = create_band(cells, cell, pixel_elev, band_size, band_map)
                    all_pixel_elevs[cell][new_band_idx].append(pixel_elev)
                    band_areas[cell][new_band_idx] += 1
                    if glacier_mask[row][col]:
                        glacier_areas[cell][new_band_idx] += 1
                        break
    # NOTE: should we check if any pixels with 0 elevation (RGM error?) are dropped 
    # into a 0-pad (invalid / nonexistent) band (since band_map will have 0 pads)?

    # Update all band median elevations for all cells, delete unused Bands
    for cell in cells:
        for band_idx, band_lower_bound in enumerate(band_map[cell]):
            cells[cell][band_idx].median_elev = np.median(all_pixel_elevs[cell][band_idx])
            # if no entries exist for this band in all_pixel_elevs, delete the band
            if not all_pixel_elevs[cell][band_idx]:
                delete_band(cells, cell, band_lower_bound, band_map)

    # Update all Band and HRU area fractions in all cells
    for cell in cells:
        for band in cell.bands:
            # Update area fraction for this Band
            new_band_area_frac = float(band_areas[cell][band_idx]) / cell_areas[cell]
            new_glacier_area_frac = float(glacier_areas[cell][band_idx]) / cell_areas[cell]
            # If the glacier HRU area fraction has changed for this band then we need to update all area fractions
            if new_glacier_area_frac != cells[cell][band_idx].area_frac_glacier:
                if new_glacier_area_frac < 0:
                    raise Exception(
                            'Calculated a negative glacier area fraction for '
                            'cell {}, band {}'.format(cell, band)
                    )
                # Calculate new non-glacier area fraction for this band
                new_non_glacier_area_frac = new_band_area_frac - new_glacier_area_frac
                # Calculate new_residual area fraction
                new_residual_area_frac = new_non_glacier_area_frac - cells[cell][band_idx].area_frac_non_glacier
# NOTE: what should be done here for Bands just created above?  cells[cell][band_idx].area_frac_non_glacier = 0 & cells[cell][band_idx].area_frac_open_ground = 1 ?
                # Calculate new open ground area fraction
                new_open_ground_area_frac = np.max([0, (cells[cell][band_idx].area_frac_open_ground + new_residual_area_frac)])
                
                # Use old proportions of vegetated areas for scaling their area fractions in this iteration
                veg_scaling_divisor = cells[cell][band_idx].area_frac_non_glacier - cells[cell][band_idx].area_frac_open_ground
                # Calculate the change in sum of vegetated area fractions
                delta_area_vegetated = np.min([0, (cells[cell][band_idx].area_frac_open_ground + new_residual_area_frac)])

                # Update area fractions for all HRUs in this Band  
                glacier_found = False
                open_ground_found = False
                for hru_idx, hru in enumerate(cells[cell][band_idx].hrus):
                    if hru.veg_type == glacier_id:
                        glacier_found = True
                        cells[cell][band_idx].hrus[hru_idx].area_frac = new_glacier_area_frac
                    if hru.veg_type == open_ground_id:
                        open_ground_found = True
                        # If open ground area fraction was reduced to 0, we delete the HRU
                        if new_open_ground_area_frac == 0:
                            band.delete_hru(open_ground_id)
                        else:
                            cells[cell][band_idx].hrus[hru_idx].area_frac = new_open_ground_area_frac
                    else:
                        # Calculate change in HRU area fraction & update
                        delta_area_hru = delta_area_vegetated * (cells[cell][band_idx].hrus[hru_idx].area_frac / veg_scaling_divisor)
                        new_hru_area_frac = cells[cell][band_id].hrus[hru_idx].area_frac + delta_area_hru
                        if new_hru_area_frac <= 0: # HRU has disappeared, never to return (only open ground can come back in its place)
                            veg_type = cells[cell][band_idx].hrus[hru_idx]
                            band.delete_hru(veg_type)
                        else:
                            cells[cell][band_idx].hrus[hru_idx].area_frac = new_hru_area_frac
                # Create glacier HRU if it didn't already exist, if we have a non-zero new_glacier_area_frac. 
                # If glacier area fraction was reduced to 0, we leave the "shadow glacier" HRU in place for VIC
                if not glacier_found and (new_glacier_area_frac > 0):
                    band.create_hru(glacier_id, new_glacier_area_frac, glacier_root_zone_parms)
                # If open ground was exposed in this band we need to add an open ground HRU
                if not open_ground_found and (new_open_ground_area_frac > 0):
                    band.create_hru(open_ground_id, new_open_ground_area_frac, open_ground_root_zone_parms)

                # Sanity check that glacier + non-glacier area fractions add up to Band's total area fraction, within tolerance
                sum_test = new_glacier_area_frac + new_non_glacier_area_frac
                if not np.allclose([sum_test], [new_band_area_frac]): # default rtol=1e-5, atol=1e-8
                    raise Exception(
                            'cell {}, band {}: glacier area fraction {} + '
                            'non-glacier area fraction {} = {} is not equal to '
                            'the band area fraction of {}'
                            .format(cell, band, new_glacier_area_frac,
                                    new_non_glacier_area_frac, sum_test,
                                    new_band_area_frac)
                    )
