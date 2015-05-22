'''vegparams.py

   This module provides a single class for reading/writing the "Vegetation
   Parameter File", a delightfully custom data serialization format, and
   provides member functions to update the vegetation parameters at runtime
   as per changes due to glacier evolution vis-a-vis RGM output.

   One can use the module as such:

   >>> vp = VegParams(glacier_id, glacier_root_zone_parms, open_ground_id, open_ground_root_zone_parms)

   Then the following attributes will be accessible:

   >>> vp.cell_ids, vp.num_veg_tiles, vp.cells
'''

from collections import OrderedDict
import csv

__all__ = ['VegParams', 'VegTileParams', 'BandInfo']

class BandInfo(object):
    """ Class capturing vegetation parameters at the elevation band level
    """
    def __init__(self, area_frac_glacier, area_frac_open_ground, area_frac_non_glacier, veg_type, area_frac, root_zone_parms ):
        self.area_frac_glacier = area_frac_glacier # this is per-band, and should have an initialization function, and gets updated during update()
        self.area_frac_non_glacier = area_frac_non_glacier
        self.area_frac_open_ground = area_frac_open_ground # this is per-band, and should have an initialization function, and gets updated during update()
        self.veg_tile_parms = [ VegTileParams(veg_type, area_frac, root_zone_parms) ]

class VegTileParams(object):
    """ Class capturing vegetation parameters at the single vegetation tile (HRU) level 
        (of which there can be many per band)
    """
    def __init__(self, veg_type, area_frac, root_zone_parms):
        self.veg_type = veg_type
        self.area_frac = area_frac
        self.root_zone_parms = root_zone_parms
                   
class VegParams(object):
    """ Class capturing all vegetation parameters information needed in the coupled VIC-RGM model
    """
    def __init__(self, veg_parm_file, glacier_id, glacier_root_zone_parms, open_ground_id, open_ground_root_zone_parms):
        self.glacier_id = glacier_id
        self.open_ground_id = open_ground_id
        self.glacier_root_zone_parms = glacier_root_zone_parms
        self.open_ground_root_zone_parms = open_ground_root_zone_parms
        self.cells = OrderedDict()
        self.area_frac_non_glacier = {}
        if veg_parm_file:
            self.load(veg_parm_file)

    @property
    def cell_ids(self):
        '''A list of all cell ids in the vegetation parameter file
        '''
        return list(self.cells.keys())

    @property
    def num_veg_tiles(self):
        '''dictionary reporting the total number vegetation tiles over each
           cell. This number is the sum of all vegetation tiles in all bands.
        '''
        count = OrderedDict()
        for cell in self.cells:
            count[cell] = 0
            for band in self.cells[cell]:
                for tile in self.cells[cell][band].veg_tile_parms:
                    count[cell]+=1
        return count

#    @staticmethod
    def read_one_cell(f):
        '''read all data from one cell and advance the file pointer to the next
        '''
        try:
            cell_id, num_veg = f.readline().split()
        except ValueError:
            return None
        cell = OrderedDict()
        for _ in range(int(num_veg)):
            line = f.readline()
            split_line = line.split()
            veg_type = split_line[0]
            area_frac = split_line[1]
            root_zone_parms = split_line[2:8]
            band_id = split_line[8]
            # If glacier_root_zone_parms or open_ground_root_zone_parms were provided at command line, 
            # override root_zone_parms here:
            if self.glacier_root_zone_parms and (veg_type == glacier_id):
                root_zone_parms = self.glacier_root_zone_parms
            elif self.open_ground_root_zone_parms and (veg_type == open_ground_id):
                root_zone_parms = self.open_ground_root_zone_parms
            try: # Create another vegetation tile for this band
                cell[band_id].veg_tile_parms.append(VegTileParams(veg_type, area_frac, root_zone_parms))
            except KeyError:
                # Create the BandInfo object for this band, along with first vegetation tile
                cell[band_id] = BandInfo(0, 0, None, veg_type, area_frac, root_zone_parms)
        return cell_id, cell

    def load(self, filename):
        '''iterate over each block of cells in the file and
           populate the class dict
        '''
        with open(filename, 'r') as f:
            while True:
                try:
                    id_, cell = self.read_one_cell(f)
                except TypeError:
                    return
                else:
                    self.cells[id_] = cell    

    def save(self, filename):
        ''' write the vegetation parameters out to a file of the same 
            format as the original vegetation parameters file
        '''
        with open(filename, 'w') as f:
            writer = csv.writer(f, delimiter=' ')
            for cell in self.cells:
                writer.writerow([cell, self.num_veg_tiles[cell]])
                for band in self.cells[cell]:
                    for tile_idx, tile in enumerate(self.cells[cell][band].veg_tile_parms):
                        line = [ self.cells[cell][band].veg_tile_parms[tile_idx].veg_type, \
                            self.cells[cell][band].veg_tile_parms[tile_idx].area_frac ] \
                            + self.cells[cell][band].veg_tile_parms[tile_idx].root_zone_parms \
                            + [ band ]
                        writer.writerow(line)

    def create_tile(self, cell_id, band_id, veg_type, area_frac, root_zone_parms):
        """ Creates a new vegetation tile of veg_type within a given cell and band in a VegParams object
        """
        new_tile = VegTileParams(veg_type, area_frac, root_zone_parms)
        # Append new_tile to existing list of tiles for this band, and sort the list ascending by veg_type
        self.cells[cell_id][band_id].veg_tile_parms.append(new_tile)
        self.cells[cell_id][band_id].veg_tile_parms.sort(key=lambda x: x.veg_type)

    def delete_tile(self, cell_id, band_id, veg_type):
        """ Deletes a vegetation tile of veg_type within a given cell and band from a VegParams object
        """
        for tile_idx, tile in enumerate(self.cells[cell_id][band_id].veg_tile_parms):
            if self.cells[cell_id][band_id].veg_tile_parms[tile_idx].veg_type == veg_type:
                del self.cells[cell_id][band_id].veg_tile_parms[tile_idx]
                break # there will only ever be one tile of any given vegetation type

    def update_tiles(cell_id, band_id, delta_area_vegetated, veg_scaling_divisor):
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
                    delete_tile(cell_id, band_id, self.open_ground_id)
                open_ground_updated = True
            # Update all other vegetation type tiles
            else:
                # Calculate change in tile area fraction & update
                delta_area_tile = delta_area_vegetated * (self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac / veg_scaling_divisor)
                new_tile_area_frac = self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac + delta_area_tile
                self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac = new_tile_area_frac
                if new_tile_area_frac <= 0: # vegetation tile has disappeared, never to return
                    delete_tile(cell_id, band_id, self.cells[cell_id][band_id].veg_tile_parms[tile].veg_type)
        # If glacier grew into this band on this iteration we need to add a glacier tile
        if not glacier_updated and (self.cells[cell_id][band_id].area_frac_glacier > 0):
            create_tile(cell_id, band_id, self.glacier_id, self.cells[cell_id][band_id].area_frac_glacier, self.glacier_root_zone_parms)
        # If open ground was exposed in this band on this iteration we need to add an open ground tile
        if not open_ground_updated and (self.cells[cell_id][band_id].area_frac_open_ground > 0):
            create_tile(cell_id, band_id, self.open_ground_id, self.cells[cell_id][band_id].area_frac_open_ground, self.open_ground_root_zone_parms)
