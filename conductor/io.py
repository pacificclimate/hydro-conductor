''' io.py 

    This module contains utility functions related to file input/output for the vic_rgm_conductor.

'''

import numpy as np

def get_rgm_pixel_mapping(pixel_map_file):
    """ Parses the RGM pixel to VIC grid cell mapping file and initialises a 2D
        grid of dimensions num_rows_dem x num_cols_dem (matching the RGM pixel
        grid), each element containing a list with the VIC cell ID associated
        with that RGM pixel and its median elevation
    """
    cell_areas = {}
    headers = {}
    with open(pixel_map_file, 'r') as f:
        # Read the number of columns and rows (order is unimportant)
        for _ in range(2):
            key, value = f.readline().split(None, 1)
            headers[key] = value
        nx = int(headers['NCOLS'])
        ny = int(headers['NROWS'])
        # create an empty two dimensional array
        cellid_map = np.empty((ny, nx))
        cellid_map.fill(np.nan)
        z_map = np.empty((ny, nx))
        z_map.fill(np.nan)
        _ = f.readline() # Consume the column headers
        for line in f:
            _, i, j, _, median_elev, cell_id = line.split()
            i, j = int(i), int(j)
            if cell_id != 'NA': #otherwise we leave it as np.NaN
                cellid_map[i,j] = cell_id
                z_map[i, j] = median_elev
            # Increment the pixel-granularity area within the grid cell
            if cell_id in cell_areas:
                cell_areas[cell_id] += 1
            else:
                cell_areas[cell_id] = 1
    return cellid_map, z_map, cell_areas, nx, ny