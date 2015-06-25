"""
    This is a set of test fixtures for the VIC-RGM Conductor.
    They are mostly based upon a simple domain consisting of a 8x8 pixel grid  
    for each VIC cell, using 3 HRU (aka vegetation) types (tree = 11, 
    open ground = 19, and glacier = 22), and a maximum of 5 elevation (aka snow) bands.  

    The initial breakdown of the first cell (ID '12345') is as follows, where
    pixels are labeled as O = open ground, T = tree, or G = glacier.  Elevation 
    bands (starting at 2000m and incrementing by a band_size of 100m) are spatially  
    comprised of concentric boxes of one pixel width, with the highest band 
    occupying the centre 4 pixels of the 8x8 grid (as open ground sticking out above 
    the glacier, represented here by lowercase 'o's).

    cell '12345':
    Spatial layout      Glacier mask

    O O O O O O O O     0 0 0 0 0 0 0 0
    O G G G G G O O     0 1 1 1 1 1 0 0 
    O G G G G G O T     0 1 1 1 1 1 0 0 
    O G G o o G O T     0 1 1 0 0 1 0 0 
    O G G o o G O T     0 1 1 0 0 1 0 0 
    O T O O O O O T     0 0 0 0 0 0 0 0 
    O T T T O O O T     0 0 0 0 0 0 0 0 
    O T T T T T T T     0 0 0 0 0 0 0 0 

    Initial HRU area fractions are calculated by adding up the sum of pixels for
    each given HRU type within a band and dividing by 64 (e.g. band 0 has a tree
    area fraction of 12/64 = 0.1875).

    The initial breakdown of the second cell (ID '23456') is as follows. Elevation 
    bands (starting at 1900m and incrementing by a band_size of 100m) are spatially 
    comprised in this domain of concentric boxes of one pixel width, with 
    the highest band / peak occupying the centre 16 (4x4) pixels of the 8x8 grid 
    (as a glacier plateau). This grid cell is located immediately to the right of 
    cell '12345' (its leftmost pixels are adjacent to the rightmost pixels of '12345').

    cell '23456':
    Spatial layout      Glacier mask

    O O G G O O O O     0 0 1 1 0 0 0 0      
    O T G G O O O O     0 0 1 1 0 0 0 0 
    T T O G G G O T     0 0 0 1 1 1 0 0 
    T T O G G G T T     0 0 0 1 1 1 0 0 
    T T O G G O T T     0 0 0 1 1 0 0 0 
    T T O O O O T T     0 0 0 0 0 0 0 0 
    T T T O O O O T     0 0 0 0 0 0 0 0 
    T T T O O T T T     0 0 0 0 0 0 0 0 

"""

import io
from pkg_resources import resource_stream, resource_filename
import numpy as np

import pytest

from conductor.cells import *
from conductor.snbparams import load_snb_parms
from conductor.vegparams import load_veg_parms

def pytest_report_header(config):
    return "VIC-RGM Conductor - Automated Test Suite"

def pytest_runtest_makereport(item, call):
    if "incremental" in item.keywords:
        if call.excinfo is not None:
            parent = item.parent
            parent._previousfailed = item

def pytest_runtest_setup(item):
    if "incremental" in item.keywords:
        previousfailed = getattr(item.parent, "_previousfailed", None)
        if previousfailed is not None:
            pytest.xfail("previous test failed (%s)" %previousfailed.name)

@pytest.fixture(scope="module")
def sample_global_file_string():
    stream = resource_stream('conductor', 'tests/input/global_peyto.txt')
    return io.TextIOWrapper(stream)

@pytest.fixture(scope="module")
def simple_unit_test_parms():

    # initial median band elevations
    test_median_elevs_simple = [2040, 2120, 2250, 2330]
    test_median_elevs = {'12345': [2040, 2120, 2250, 2330],
                        '23456': [1970, 2005, 2120]} 

    # Just for single cell unit tests:
    test_area_fracs_simple = [0.1875, 0.25, # Band 0 (11, 19)
                        0.0625, 0.125, 0.125, # Band 1 (11, 19, 22)
                        0.0625, 0.125, # Band 2 (19, 22)
                        0.0625] # Band 3 (19)
    test_area_fracs = {'12345': [0.1875, 0.25, # Band 0 (11, 19)
                        0.0625, 0.125, 0.125, # Band 1 (11, 19, 22)
                        0.0625, 0.125, # Band 2 (19, 22)
                        0.0625], # Band 3 (19)
                    '23456': [0.25, 0.15625, 0.03125, # Band 1 (11, 19, 22)
                        0.15625, 0.125, 0.03125, # Band 2 (11, 19, 22)
                        0.125, 0.125]} # Band 3 (19, 22)
    test_area_fracs_by_band = {'12345': {'0': [0.1875, 0.25], # Band 0 (11, 19)
                        '1': [0.0625, 0.125, 0.125], # Band 1 (11, 19, 22)
                        '2': [0.0625, 0.125], # Band 2 (19, 22)
                        '3': [0.0625]}, # Band 3 (19)
                    '23456': {'1': [0.25, 0.15625, 0.03125], # Band 1 (11, 19, 22)
                        '2': [0.15625, 0.125, 0.03125], # Band 2 (11, 19, 22)
                        '3': [0.125, 0.125]} } # Band 3 (19, 22)

    test_veg_types = [11, 19, 22]

    expected_num_hrus = {'12345': [2, 3, 2, 1],
                '23456': [3, 3, 2] }

    expected_root_zone_parms = {'11': [0.10, 0.60, 0.20, 0.25, 1.70, 0.15], # 11
                        '19': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0], # 19
                        '22': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0]} # 22

    return test_median_elevs_simple, test_median_elevs, test_area_fracs_simple, test_area_fracs, \
            test_area_fracs_by_band, test_veg_types, expected_num_hrus, expected_root_zone_parms

@pytest.fixture(scope="module")
def large_merge_cells_unit_test_parms():
    fname = resource_filename('conductor', 'tests/input/snow_band.txt')
    elevation_cells = load_snb_parms(fname, 15)
    fname = resource_filename('conductor', 'tests/input/veg.txt')
    hru_cells = load_veg_parms(fname)
    expected_zs = [ 2076, 2159, 2264, 2354, 2451, 2550, 2620, 2714, 2802 ]
    expected_afs = {0.000765462339, 0.000873527611, 0.009125511809, 0.009314626034, 0.004426673711, 0.004558753487, 0.001388838859, 0.000737445417}

    return elevation_cells, hru_cells, expected_zs, expected_afs

@pytest.fixture(scope="module")
def toy_domain_64px_cells():

    # NOT USED BUT USEFUL INFO: initial band map of lower elevation bounds for all existing (valid) bands
    # test_band_map = {'12345': [2000, 2100, 2200, 2300, 0], # allows for glacier growth at top
    #            '23456': [0, 1900, 2000, 2100, 0]} # allows for glacier growth at top, and revelation of lower band at bottom

    fname = resource_filename('conductor', 'tests/input/snb_toy_64px.txt')
    elevation_cells = load_snb_parms(fname, 5)
    fname = resource_filename('conductor', 'tests/input/vpf_toy_64px.txt')
    hru_cells = load_veg_parms(fname)
    cells = merge_cell_input(hru_cells, elevation_cells)
    cell_ids = list(cells.keys())
    # We have a total allowable number of snow bands of 5, with 100m spacing
    num_snow_bands = 5
    band_size = 100
    # Initially we have just 4 bands loaded for cell 0, and 3 for cell 1
    expected_band_ids = {cell_ids[0]: [0, 1, 2, 3],
                        cell_ids[1]: [1, 2, 3]}
    expected_root_zone_parms = {'11': [0.10, 0.60, 0.20, 0.25, 1.70, 0.15], # 11
                        '19': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0], # 19
                        '22': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0]} # 22

    # Spatial DEM layout
    initial_dem_by_cells = { cell_ids[0]:                  
                                np.array([[2065, 2055, 2045, 2035, 2025, 2015, 2005, 2000],
                                        [2075, 2100, 2120, 2140, 2130, 2120, 2100, 2005],
                                        [2085, 2110, 2250, 2270, 2260, 2240, 2110, 2010],
                                        [2090, 2120, 2260, 2377, 2310, 2245, 2125, 2015],
                                        [2070, 2110, 2250, 2340, 2320, 2250, 2130, 2020],
                                        [2090, 2105, 2200, 2210, 2220, 2220, 2120, 2015],
                                        [2090, 2100, 2105, 2110, 2140, 2150, 2130, 2010],
                                        [2080, 2075, 2065, 2055, 2045, 2035, 2020, 2000] ]),

                            cell_ids[1]:
                                np.array([[1970, 1975, 1995, 1995, 1975, 1965, 1960, 1960],
                                        [1970, 2000, 2045, 2055, 2005, 2005, 2000, 1965],
                                        [1975, 2000, 2100, 2155, 2160, 2140, 2000, 1970],
                                        [1985, 2005, 2105, 2160, 2180, 2130, 2000, 1975],
                                        [1990, 2010, 2110, 2150, 2140, 2105, 2005, 1980],
                                        [1980, 2005, 2105, 2105, 2110, 2100, 2000, 1980],
                                        [1970, 2000, 2000, 2020, 2035, 2025, 2000, 1970],
                                        [1965, 1965, 1970, 1970, 1975, 1960, 1950, 1960] ])
                    }

    initial_glacier_mask_by_cells = { cell_ids[0]: 
                                        np.array([[ 0, 0, 0, 0, 0, 0, 0, 0 ],
                                                  [ 0, 1, 1, 1, 1, 1, 0, 0 ],
                                                  [ 0, 1, 1, 1, 1, 1, 0, 0 ],
                                                  [ 0, 1, 1, 0, 0, 1, 0, 0 ],
                                                  [ 0, 1, 1, 0, 0, 1, 0, 0 ],
                                                  [ 0, 0, 0, 0, 0, 0, 0, 0 ],
                                                  [ 0, 0, 0, 0, 0, 0, 0, 0 ],
                                                  [ 0, 0, 0, 0, 0, 0, 0, 0 ]]),
                                    cell_ids[1]:
                                        np.array([[ 0, 0, 1, 1, 0, 0, 0, 0 ],
                                                  [ 0, 0, 1, 1, 0, 0, 0, 0 ],
                                                  [ 0, 0, 0, 1, 1, 1, 0, 0 ],
                                                  [ 0, 0, 0, 1, 1, 1, 0, 0 ],
                                                  [ 0, 0, 0, 1, 1, 0, 0, 0 ],
                                                  [ 0, 0, 0, 0, 0, 0, 0, 0 ],
                                                  [ 0, 0, 0, 0, 0, 0, 0, 0 ],
                                                  [ 0, 0, 0, 0, 0, 0, 0, 0 ]])
                                    }

    def build_padded_dem_aligned_map(array_1, array_2, padding_thickness):
        ''' Helper function to create a map that is aligned with a DEM made up of two adjacent
            rectangular pixel arrays (representing 2 VIC cells), with NaN padding around
            the edges (simulating the output from get_rgm_pixel_mapping(), where some
            pixels that are read in do not belong to either VIC cell). This can be used to 
            generate test instances of surf_dem, bed_dem, cellid_map, glacier_mask
        '''
        vertical_pad = np.empty((padding_thickness, 2*padding_thickness+len(array_1[0])+len(array_2[0])))
        vertical_pad.fill(np.nan)
        horizontal_pad = np.empty((len(array_1), padding_thickness))
        horizontal_pad.fill(np.nan)
        padded_map = np.concatenate((array_1, array_2), axis=1)
        padded_map = np.concatenate((horizontal_pad, padded_map), axis=1)
        padded_map = np.concatenate((padded_map, horizontal_pad), axis=1)
        padded_map = np.concatenate((vertical_pad, padded_map), axis=0)
        padded_map = np.concatenate((padded_map, vertical_pad), axis=0)
        return padded_map

    # Create cellid_map that is padded by 2 rows & columns of np.NaNs (to simulate mapping
    # of valid pixels to VIC cells, as read in by get_rgm_pixel_mapping())
    cellid_map_cell_0 = np.empty((8,8))
    cellid_map_cell_0.fill(cell_ids[0])
    cellid_map_cell_1 = np.empty((8,8))
    cellid_map_cell_1.fill(cell_ids[1])
    cellid_map = build_padded_dem_aligned_map(cellid_map_cell_0, cellid_map_cell_1, 2)

    # Create initial surf_dem with padding of 2
    surf_dem = build_padded_dem_aligned_map(initial_dem_by_cells[cell_ids[0]], initial_dem_by_cells[cell_ids[1]], 2)
    # Non-padded version of surf_dem:
    #surf_dem = np.concatenate((initial_dem_by_cells[cell_ids[0]], initial_dem_by_cells[cell_ids[1]]), axis=1)

    # Create initial glacier mask with padding of 2
    glacier_mask = build_padded_dem_aligned_map(initial_glacier_mask_by_cells[cell_ids[0]], initial_glacier_mask_by_cells[cell_ids[1]], 2)

    # The toy DEM above, broken down by elevation bands.  Useful for checking median elevations
    cell_band_pixel_elevations = {
        cell_ids[0]: [
                    # band 0, median: 2040
                    [2065, 2055, 2045, 2035, 2025, 2015, 2005, 2000,
                     2075,                                     2005,
                     2085,                                     2010,
                     2090,                                     2015,
                     2070,                                     2020,
                     2090,                                     2015,
                     2090,                                     2010,
                     2080, 2075, 2065, 2055, 2045, 2035, 2020, 2000],

                    # band 1, median: 2120
                    [2100, 2120, 2140, 2130, 2120, 2100,
                     2110,                         2110,
                     2120,                         2125,
                     2110,                         2130,
                     2105,                         2120,
                     2100, 2105, 2110, 2140, 2150, 2130],

                    # band 2, median: 2250
                    [2250, 2270, 2260, 2240,
                     2260,             2250,
                     2250,             2250,
                     2200, 2210, 2220, 2220],

                    # band 3, median: 2330
                    [2377, 2310,
                    2340, 2320]
                    ],
        cell_ids[1]: [
                    # band 0, median: 1970
                    [1970, 1975, 1995, 1995, 1975, 1965, 1960, 1960,
                    1970,                                       1965,
                    1975,                                       1970,
                    1985,                                       1975,
                    1990,                                       1980,
                    1980,                                       1980,
                    1970,                                       1970,
                    1965, 1965, 1970, 1970, 1975, 1960, 1950, 1960],

                    # band 1, median: 2005
                    [2000, 2045, 2055, 2005, 2005, 2000,
                    2000,                           2000,
                    2005,                           2000,
                    2010,                           2005,
                    2005,                           2000,
                    2000, 2000, 2020, 2035, 2025, 2000],

                    # band 2, median: 2120
                    [2100, 2155, 2160, 2140,
                    2105, 2160, 2180, 2130,
                    2110, 2150, 2140, 2105,
                    2105, 2105, 2110, 2100],
                    ] }

    return cells, cell_ids, num_snow_bands, band_size, expected_band_ids, expected_root_zone_parms, \
            cellid_map, surf_dem, glacier_mask, cell_band_pixel_elevations

