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

__all__ = ['HydroResponseUnit', 'Band']


def load_veg_parms(veg_parm_file, glacier_id, open_ground_id, glacier_root_zone_parms=None, open_ground_root_zone_parms=None):
    """ Reads in VIC vegetation parameter file and creates and partially 
        initializes all VIC grid cells.
    """
    def load(filename):
        '''iterate over each block of cells in the file and
           populate the class dict
        '''
        cells = OrderedDict()
        with open(filename, 'r') as f:
            while True:
                try:
                    id_, cell = read_one_cell(f)
                except TypeError:
                    break
                else:
                    cells[id_] = cell
        return cells
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
            veg_type = int(split_line[0])
            area_frac = float(split_line[1])
            root_zone_parms = [float(x) for x in split_line[2:8]]
            band_id = split_line[8]
            # If glacier_root_zone_parms or open_ground_root_zone_parms were provided at command line, 
            # override root_zone_parms here:
            if (veg_type == glacier_id) and glacier_root_zone_parms is not None:
                root_zone_parms = glacier_root_zone_parms
            elif (veg_type == open_ground_id) and open_ground_root_zone_parms is not None:
                root_zone_parms = open_ground_root_zone_parms
            try: # Create another HRU for this band
                cell[band_id].hrus.append(HydroResponseUnit(veg_type, area_frac, root_zone_parms))
            except KeyError:
                # Create the BandInfo object for this band, along with first vegetation tile
                cell[band_id] = Band(area_frac, None, glacier_id, open_ground_id)
                cell[band_id].hrus.append(HydroResponseUnit(veg_type, area_frac, root_zone_parms))
        return cell_id, cell
    
    cells = load(veg_parm_file)
    return cells

def save_veg_parms(cells, filename):
    ''' write the vegetation parameters out to a file of the same 
        format as the original vegetation parameters file
    '''
    with open(filename, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        for cell in cells:
            writer.writerow([cell, sum([cells[cell][band].num_hrus for band in cells[cell]])])
            for band in cells[cell]:
                for hru in cells[cell][band].hrus:
                    line = [hru.veg_type, hru.area_frac] + hru.root_zone_parms + [band]
                    writer.writerow(line)

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

class HydroResponseUnit(object):
    """ Class capturing vegetation parameters at the single vegetation tile (HRU) level 
        (of which there can be many per band)
    """
    def __init__(self, veg_type, area_frac, root_zone_parms):
        self.veg_type = veg_type
        self.area_frac = area_frac
        self.root_zone_parms = root_zone_parms


#TODO:  
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

def delete_band(cell_id, band_lower_bound):
    """ Removes the band starting at band_lower_bound from the band_map and 
        sets area_fracs and median_elevs to 0 pads (VIC needs the number 
        of band placeholders to remain constant)
    """
    band_idx = self.cells[cell_id].band_map.index(band_lower_bound)
    self.cells[cell_id].area_fracs[band_idx] = 0
    self.cells[cell_id].median_elevs[band_idx] = 0
    del self.cells[cell_id].band_map[band_idx]



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

                   
