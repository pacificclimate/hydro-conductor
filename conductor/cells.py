''' cells.py

    This module captures the conceptual components of the VIC grid cell,
    broken into classes by elevation band and Hydrological Response Unit (HRU) 

'''

class Band(object):
    """ Class capturing vegetation parameters at the elevation band level
    """
    def __init__(self, area_frac, median_elev, glacier_id, open_ground_id):
        self.median_elev = median_elev
        self.hrus = []
        self.glacier_id = glacier_id
        self.open_ground_id = open_ground_id
    @property
    def num_hrus(self):
        return len(self.hrus)
    @property
    def area_frac(self): 
        """ The Band area fraction, equal to the sum of HRU area fractions within this band, 
            which should be equal to the total area fraction of the band as given in the 
            Snow Band Parameters file
        """
        return sum([hru.area_frac for hru in self.hrus])
    @property
    def area_frac_glacier(self):
        #return sum([hru.area_frac if hru.veg_type == glacier_id for hru in self.hrus])
        return sum([hru.area_frac for hru in self.hrus if hru.veg_type == self.glacier_id])
        #return self.hrus[self.hrus.veg_type == glacier_id].area_frac
    @property
    def area_frac_non_glacier(self):
        return self.area_frac - self.area_frac_glacier

    # TODO: modify for new cell/Band implementation
    def create_band(cell_id, elevation):
        """ Creates a new elevation band for a cell with an initial median elevation
        """
        band_lower_bound = int(elevation - elevation % self.band_size)
        bisect.insort_left(self.cells[cell_id].band_map, band_lower_bound)
        band_idx = self.cells[cell_id].band_map.index(band_lower_bound)
        # area_fracs for this band will already be represented by a pad of 0.0, and will get updated in update_band_area_fracs()
        # Replace a 0 pad value with an initial median elevation, or if no spaces are left raise an error
        try:
            self.cells[cell_id].median_elevs[band_idx] = elevation  
        except IndexError:
            print('SnbParams::create_band: IndexError: Attempted to create a new elevation band \
                at {}m (RGM output DEM pixel elevation at {}) in cell {}, but ran out of available slots. \
                Increase 0 padding in the VIC Snow Band Parameters file and re-run. Exiting.\
                \n'.format(band_lower_bound, elevation, cell_id))
            sys.exit(0)
        return band_idx

    # TODO: modify for new cell/Band implementation
    def delete_band(cell_id, band_lower_bound):
        """ Removes the band starting at band_lower_bound from the band_map and 
            sets area_fracs and median_elevs to 0 pads (VIC needs the number 
            of band placeholders to remain constant)
        """
        band_idx = self.cells[cell_id].band_map.index(band_lower_bound)
        self.cells[cell_id].area_fracs[band_idx] = 0
        self.cells[cell_id].median_elevs[band_idx] = 0
        del self.cells[cell_id].band_map[band_idx]


class HydroResponseUnit(object):
    """ Class capturing vegetation parameters at the single vegetation tile (HRU) level 
        (of which there can be many per band)
    """
    def __init__(self, veg_type, area_frac, root_zone_parms):
        self.veg_type = veg_type
        self.area_frac = area_frac
        self.root_zone_parms = root_zone_parms

    #TODO: remove references to band_id, as the HRU is now nested within a band
    def create_hru(self, cell_id, band_id, veg_type, area_frac, root_zone_parms):
        """ Creates a new vegetation tile of veg_type within a given cell and band in a VegParams object
        """
        new_tile = HydroResponseUnit(veg_type, area_frac, root_zone_parms)
        # Append new_tile to existing list of tiles for this band, and sort the list ascending by veg_type
        self.cells[cell_id][band_id].veg_tile_parms.append(new_tile)
        self.cells[cell_id][band_id].veg_tile_parms.sort(key=lambda x: x.veg_type)

    #TODO: remove references to band_id, as the HRU is now nested within a band
    def delete_hru(self, cell_id, band_id, veg_type):
        """ Deletes a vegetation tile of veg_type within a given cell and band from a VegParams object
        """
        for tile_idx, tile in enumerate(self.cells[cell_id][band_id].veg_tile_parms):
            if self.cells[cell_id][band_id].veg_tile_parms[tile_idx].veg_type == veg_type:
                del self.cells[cell_id][band_id].veg_tile_parms[tile_idx]
                break # there will only ever be one tile of any given vegetation type

    #TODO: remove references to band_id, as the HRU is now nested within a band
    def update_hrus(cell_id, band_id, delta_area_vegetated, veg_scaling_divisor):
        """ Updates all vegetation tile area fractions within a band in a cell_id, 
            and creates / deletes tiles if necessary
        """
        glacier_updated = False
        open_ground_updated = False
        # Iterate through all existing vegetation tiles in this band and update their area fractions
        for tile in self.cells[cell_id][band_id].veg_tile_parms:
            # Update glacier tile, if it exists
            if self.cells[cell_id][band_id].veg_tile_parms[tile].veg_type == self.glacier_id:
                # NOTE: we allow glacier area fractions of 0 for the purposes of VIC's shadow glacier requirements.
                # However, by never deleting a glacier tile we may end up with shadow glaciers across many bands
                # by the end of a run (even if glacier only appears in some bands "momentarily") 
                self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac = self.cells[cell_id][band_id].area_frac_glacier
                glacier_updated = True
            # Update open ground tile, if it exists
            elif self.cells[cell_id][band_id].veg_tile_parms[tile].veg_type == self.open_ground_id:
                if self.cells[cell_id][band_id].area_frac_open_ground > 0:
                    self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac = self.cells[cell_id][band_id].area_frac_open_ground
                else: # remove the tile so it is not written to the veg parameters file
                    delete_hru(cell_id, band_id, self.open_ground_id)
                open_ground_updated = True
            # Update all other vegetation type tiles
            else:
                # Calculate change in tile area fraction & update
                delta_area_tile = delta_area_vegetated * (self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac / veg_scaling_divisor)
                new_tile_area_frac = self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac + delta_area_tile
                self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac = new_tile_area_frac
                if new_tile_area_frac <= 0: # vegetation tile has disappeared, never to return
                    delete_hru(cell_id, band_id, self.cells[cell_id][band_id].veg_tile_parms[tile].veg_type)
        # If glacier grew into this band on this iteration we need to add a glacier tile
        if not glacier_updated and (self.cells[cell_id][band_id].area_frac_glacier > 0):
            create_hru(cell_id, band_id, self.glacier_id, self.cells[cell_id][band_id].area_frac_glacier, self.glacier_root_zone_parms)
        # If open ground was exposed in this band on this iteration we need to add an open ground tile
        if not open_ground_updated and (self.cells[cell_id][band_id].area_frac_open_ground > 0):
            create_hru(cell_id, band_id, self.open_ground_id, self.cells[cell_id][band_id].area_frac_open_ground, self.open_ground_root_zone_parms)

                       
