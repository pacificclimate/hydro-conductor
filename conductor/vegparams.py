'''vegparams.py

   This module provides a single class for reading/writing the "Vegetation
   Parameter File", a delightfully custom data serialization format.
   One can use the module as such:

   >>> vp = VegParams(my_file)

   Then the following attributes will be accessible:

   >>> vp.cell_ids, vp.num_veg_tiles, vp.cells
'''

from collections import OrderedDict
import csv

__all__ = ['VegParams', 'VegTileParms', 'BandInfo']

class BandInfo(object):
    def __init__(self, area_frac_glacier, area_frac_open_ground, area_frac_non_glacier, veg_type, area_frac, root_zone_parms ):
        self.area_frac_glacier = area_frac_glacier # this is per-band, and should have an initialization function, and gets updated during update()
        self.area_frac_non_glacier = area_frac_non_glacier
        self.area_frac_open_ground = area_frac_open_ground # this is per-band, and should have an initialization function, and gets updated during update()
        self.veg_tile_parms = [ VegTileParms(veg_type, area_frac, root_zone_parms) ]

class VegTileParms(object):
    def __init__(self, veg_type, area_frac, root_zone_parms):
        self.veg_type = veg_type
        self.area_frac = area_frac
        self.root_zone_parms = root_zone_parms
                   
class VegParams(object):
    # NOTE: these 4 parameters should probably be set when a class object is instantiated, 
    # as Markus has specified that open_ground_id be a global parameter and the glacier/open 
    # ground root_zone_parms be options at invocation
    glacier_id = '22'
    open_ground_id ='19'
    glacier_root_zone_parms = ''
    open_ground_root_zone_parms = ''

    def __init__(self, veg_parm_file=None):
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

    @staticmethod
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
            # TODO: if glacier_root_zone_parms or open_ground_zone_parms were provided at command line, 
            # override root_zone_parms here:
            #if veg_type is glacier_id and glacier_root_zone_parms:
            #    root_zone_parms = glacier_root_zone_parms
            #elif veg_type is open_ground and open_ground_zone_parms:
            #    root_zone_parms = open_ground_zone_parms
            try:
                cell[band_id].veg_tile_parms.append(VegTileParms(veg_type, area_frac, root_zone_parms))
            except KeyError:
                # TODO: check that it is OK to set these first two init parms to 0
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
        pass

    def delete_tile(self, cell_id, band_id, veg_type):
        """ Deletes a vegetation tile of veg_type within a given cell and band from a VegParams object
        """
        pass

    def update_tiles(cell_id, band_id, delta_area_vegetated, veg_scaling_divisor):
        """ Updates all vegetation tile area fractions within a band in a cell_id, 
            and creates / deletes tiles if necessary
        """
        glacier_updated = False
        open_ground_updated = False
        for tile in self.cells[cell_id][band_id].veg_tile_parms:
            if self.cells[cell_id][band_id].veg_tile_parms[tile].veg_type == self.glacier_id:
                # we allow glacier area fractions of 0 for the purposes of VIC's shadow glacier requirements
                # NOTE: by never deleting a glacier tile we may end up with shadow glaciers across many bands
                # by the end of a run (even if glacier only appears in some bands "momentarily") 
                self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac = self.cells[cell_id][band_id].area_frac_glacier
                glacier_updated = True
            elif self.cells[cell_id][band_id].veg_tile_parms[tile].veg_type == self.open_ground_id:
                if self.cells[cell_id][band_id].area_frac_open_ground > 0:
                    self.cells[cell_id][band_id].veg_tile_parms[tile].area_frac = self.cells[cell_id][band_id].area_frac_open_ground
                        # if this fails, create an open ground tile
                else: # remove the tile so it is not written to the veg parameters file
                    delete_tile(cell_id, band_id, self.open_ground_id)
                open_ground_updated = True
            else:
        #TODO: apply delta area calculation & update for all other vegetation types


        # If glacier grew into this band on this iteration we need to add a glacier tile
        if not glacier_updated and (self.cells[cell_id][band_id].area_frac_glacier > 0):
            create_tile(cell_id, band_id, self.glacier_id, self.cells[cell_id][band_id].area_frac_glacier, self.glacier_root_zone_parms)
        # If open ground was exposed in this band on this iteration we need to add an open ground tile
        if not open_ground_updated and (self.cells[cell_id][band_id].area_frac_open_ground > 0):
            create_tile(cell_id, band_id, self.open_ground_id, self.cells[cell_id][band_id].area_frac_open_ground, self.open_ground_root_zone_parms)
