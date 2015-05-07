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

__all__ = ['VegParams']

def leaves(d):
    '''Return the number of leaves in a multilevel dict
    '''
    if hasattr(d, 'values'):
        return sum([leaves(value) for value in d.values()])
    elif hasattr(d, '__len__'):
        return len(d)
    else:
        return 1
                    
class VegParams(object):
    def __init__(self, veg_parm_file=None):
        self.cells = OrderedDict()
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
        # Count is two levels deep
        return OrderedDict( (cell_id, leaves(cell)) for cell_id, cell in self.cells.items() )

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
            band_id = line.split()[-1]
            try:
                cell[band_id].append(line.split())
            except KeyError:
                cell[band_id] = [ (line.split()) ]
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
                    for line in self.cells[cell][band]:
                        writer.writerow(line)
                        print(' '.join(map(str, line)))
