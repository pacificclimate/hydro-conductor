'''vegparams.py

   This module provides functions for reading/writing the VIC Vegetation
   Parameter File, a delightfully custom data serialization format.

'''

from collections import OrderedDict
import csv
from conductor.cells import Band, HydroResponseUnit

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
        #for cell in cells.values():
        for cell in cells:
            writer.writerow([cell, sum([cells[cell][band].num_hrus for band in cells[cell]])])
            for band in cells[cell]:
                for hru in cells[cell][band].hrus:
                    line = [hru.veg_type, hru.area_frac] + hru.root_zone_parms + [band]
                    writer.writerow(line)