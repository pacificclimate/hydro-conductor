'''snbparams.py

   This module represents the VIC "Snow Band File",
   which captures all the elevation band information
   for all VIC grid cells in a given model run.
   The format of the snow band file is one line per VIC cell:
   cell_id_0 area_frac_band_0 ... area_frac_band_N median_elev_band_0 ... median_elev_band_N
   where N should be equal to num_snow_bands
'''

__all__ = ['load_snb_parms', 'save_snb_parms']

import bisect
from collections import OrderedDict
import csv

def load_snb_parms(cells, snb_file, num_snow_bands, band_size):
    """ Reads in a Snow Band Parameter File and populates the median elevation
        attribute for an existing set of VIC cells and creates a band map to
        keep track of the lower bounds of each band and zero pads provided by 
        the user in the file, to allow for glacier growth/slide into 
        previously non-existent elevations (VIC requires the zero pads).
    """
    with open(snb_file, 'r') as f:
        median_elevs = OrderedDict()
        band_map = OrderedDict()
        for line in f:
            #print('snb file line: {}'.format(line))
            split_line = line.split()
            cell_id = split_line[0]
            num_columns = len(split_line)
            # Should have the cell_id followed by num_snow_bands columns 
            # for each of area fractions, median elevations, and optionally, Pfactors:
            if num_columns % num_snow_bands != 1:
                print('SnbParams::load(): Error: Number of columns ({}) in snow band file {} \
                    is incorrect for the given number of snow bands ({}) given in the global \
                    parameter file (should be a multiple of num_snow_bands, plus 1). \
                    Exiting.\n'.format(num_columns, snb_file, num_snow_bands))
                sys.exit(0)
            median_elevs[cell_id] = [int(x) for x in split_line[num_snow_bands+1 : 2*num_snow_bands+1]]
    # Create and populate the band_map for keeping track of existing and zero-pad bands,
    # and set the median_elevation property for each existing band
    for cell in cells:
        band_map[cell] = []
        for band_idx, elev in enumerate(median_elevs[cell]):
            band_map[cell].append(int(elev - elev % band_size))
            if elev > 0: # there is no Band object within the cell for zero pad elevations
                cells[cell][str(band_idx)].median_elev = elev
    return band_map

def save_snb_parms(cells, filename, band_map):
    """ Writes current (updated) snow band parameters to a new temporary
        Snow Band File for feeding back into VIC
    """
    with open(filename, 'w') as f:
         writer = csv.writer(f, delimiter=' ')
         for cell in cells:
            line = [cell]
            for band_idx, band in enumerate(band_map[cell]):
                if band_map[cell][band_idx] != 0:
                    line.append(cells[cell][str(band_idx)].area_frac)
                else:
                    line.append(0)  
            for band_idx, band in enumerate(band_map[cell]):
                if band_map[cell][band_idx] != 0:
                    line.append(cells[cell][str(band_idx)].median_elev)
                else:
                    line.append(0)
            writer.writerow(line)
