'''snbparams.py

   This module represents the VIC "Snow Band File",
   which captures all the elevation band information
   for all VIC grid cells in a given model run.
   The format of the snow band file is one line per VIC cell:
   cell_id_0 area_frac_band_0 ... area_frac_band_N median_elev_band_0 ... median_elev_band_N Pfactor_band_0 ... Pfactor_band_N
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
        '''A list of all cell ids in the vegetation parameter file
        '''
        return list(self.cells.keys())

    def create_cell(self, cell_id, area_fracs, median_elevs, pfactors):
        self.cells[cell_id] = namedtuple('area_fracs', 'median_elevs', 'pfactors', 'band_map')
        self.cells[cell_id].area_fracs = area_fracs
        self.cells[cell_id].median_elevs = median_elevs
        self.cells[cell_id].pfactors = pfactors
        self.cells[cell_id].band_map = []
        # create a list of elevation bands for the grid cell of band_size width in meters, rounded down to nearest increment of band_size
        for elev in median_elevs:
            self.cells[cell_id].band_map.append(int(elev - elev % self.band_size))

    def create_band(self, cell_id, pixel_elev):
        '''
        Creates a new elevation band in the band_map for cell_id, rounded down to the nearest increment of band_size
        '''
        bisect.insort_left(self.cells[cell_id].band_map, int(pixel_elev - pixel_elev % self.band_size))
        
    def delete_band(self, cell_id, elevation):
        '''
        Removes the band starting at elevation from the band_map for cell_id
        '''
        del self.cells[cell_id].band_map[bisect.bisect_left(self.cells[cell_id].band_map, elevation)]

    def update_area_frac(self, band):
        pass

    def load(self, snb_file):
        ''' Reads in a Snow Band File
        '''
        with open(snb_file, 'r') as f:
            for line in f:
                #print('snb file line: {}'.format(line))
                split_line = line.split()
                num_columns = len(split_line)
                cell_id = split_line[0]
                if num_columns != 3*self.num_snow_bands + 1:
                    print('SnbParams::load(): Error: Number of columns ({}) in snow band file {} is incorrect for the given number of snow bands ({}) given in the global parameter file (should be 3 * num_snow_bands + 1). Exiting.\n'.format(num_columns, snb_file, self.num_snow_bands))
                    sys.exit(0)
                self.create_cell(cell_id, [float(x) for x in split_line[1 : self.num_snow_bands+1]],[int(x) for x in split_line[self.num_snow_bands+1 : 2*self.num_snow_bands+1]],[float(x) for x in split_line[2*self.num_snow_bands+1 : 3*self.num_snow_bands+1]])

    def save(self, filename):
        """ Writes current (updated) snow band parameters to a new temporary
            Snow Band File for feeding back into VIC
        """
        with open(filename, 'w') as f:
             writer = csv.writer(f, delimiter=' ')
             for cell in self.cells:                  
                 line = [cell] + self.cells[cell].area_fracs + self.cells[cell].median_elevs + self.cells[cell].pfactors
                 writer.writerow(line)
