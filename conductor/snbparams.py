'''snbparams.py

   This module represents the VIC "Snow Band File",
   which captures all the elevation band information
   for all VIC grid cells in a given model run.
   The format of the snow band file is one line per VIC cell:
   cell_id_0 area_frac_band_0 ... area_frac_band_N median_elev_band_0 ... median_elev_band_N
   where N should be equal to num_snow_bands
'''

__all__ = ['get_snb_parms', 'create_band_map', 'update_snb_parms', 'write_snb_parms_file']

import bisect
from collections import namedtuple, OrderedDict
import csv
import numpy as np

class SnbParams(object):

    def __init__(self, num_snow_bands, band_size, snb_parm_file=None):
        self.num_snow_bands = num_snow_bands
        self.band_size = band_size
        self.cells = OrderedDict()
        if snb_parm_file:
            self.load(snb_parm_file)

    @property
    def cell_ids(self):
        """ Returns a list of all cell ids in the vegetation parameter file
        """
        return list(self.cells.keys())

    def create_cell(self, cell_id, area_fracs, median_elevs):
        """ Creates an entry for a VIC grid cell from a set of area fractions, and median elevations
        """
        self.cells[cell_id] = namedtuple('area_fracs', 'median_elevs', 'band_map')
        self.cells[cell_id].area_fracs = area_fracs
        self.cells[cell_id].median_elevs = median_elevs
        self.cells[cell_id].band_map = [] # list of valid elevation bands
        # create a list of elevation bands for the grid cell of band_size width in meters, rounded down to nearest increment of band_size
        for elev in median_elevs:
            if elev > 0: # we ignore any 0 padding representing possible bands that don't (yet) exist
                self.cells[cell_id].band_map.append(int(elev - elev % self.band_size))

    def create_band(self, cell_id, elevation):
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

    def delete_band(self, cell_id, band_lower_bound):
        """ Removes the band starting at band_lower_bound from the band_map and 
            sets area_fracs and median_elevs to 0 pads (VIC needs the number 
            of band placeholders to remain constant)
        """
        band_idx = self.cells[cell_id].band_map.index(band_lower_bound)
        self.cells[cell_id].area_fracs[band_idx] = 0
        self.cells[cell_id].median_elevs[band_idx] = 0
        del self.cells[cell_id].band_map[band_idx]

    def load(self, snb_file):
        """ Reads in a Snow Band Parameter File
        """
        with open(snb_file, 'r') as f:
            for line in f:
                #print('snb file line: {}'.format(line))
                split_line = line.split()
                num_columns = len(split_line)
                cell_id = split_line[0]
                if num_columns != 3*self.num_snow_bands + 1:
                    print('SnbParams::load(): Error: Number of columns ({}) in snow band file {} is incorrect for the given number of snow bands ({}) given in the global parameter file (should be 3 * num_snow_bands + 1). Exiting.\n'.format(num_columns, snb_file, self.num_snow_bands))
                    sys.exit(0)
                self.create_cell(cell_id, [float(x) for x in split_line[1 : self.num_snow_bands+1]],[int(x) for x in split_line[self.num_snow_bands+1 : 2*self.num_snow_bands+1]]])

    def save(self, filename):
        """ Writes current (updated) snow band parameters to a new temporary
            Snow Band File for feeding back into VIC
        """
        with open(filename, 'w') as f:
             writer = csv.writer(f, delimiter=' ')
             for cell in self.cells:
#                if len(self.cells[cell].area_fracs) < self.num_snow_bands:

                line = [cell] + self.cells[cell].area_fracs + self.cells[cell].median_elevs + self.cells[cell].pfactors
                writer.writerow(line)
