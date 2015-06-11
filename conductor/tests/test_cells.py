''' This is a set of tests for the cells.py module.
    It is based upon a simple 8x8 pixel grid domain per VIC cell, 
    using 3 HRU (aka vegetation) types (tree = 11, open ground = 19, 
    and glacier = 22), and a maximum of 5 elevation (aka snow) bands.  

    The initial breakdown of the first cell (ID '12345') is as follows, where
    pixels are labeled as O = open ground, T = tree, or G = glacier.  Elevation 
    bands (starting at 2000m and incrementing by a band_size of 100m) are spatially  
    comprised in this domain of concentric boxes of one pixel width, with 
    the highest band / peak occupying the centre 4 pixels of the 8x8 grid 
    (as open ground sticking out above the glacier, shown by lowercase 'o's).

    cell '12345'
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

    O O G G O O O O     0 0 1 1 0 0 0 0      
    O T G G O O O O     0 0 1 1 0 0 0 0 
    T T O G G G O T     0 0 0 1 1 1 0 0 
    T T O G G G T T     0 0 0 1 1 1 0 0 
    T T O G G O T T     0 0 0 1 1 0 0 0 
    T T O O O O T T     0 0 0 0 0 0 0 0 
    T T T O O O O T     0 0 0 0 0 0 0 0 
    T T T O O T T T     0 0 0 0 0 0 0 0 

'''

from collections import OrderedDict
from pkg_resources import resource_filename

import pytest

from conductor.cells import *
from conductor.snbparams import load_snb_parms, PaddedDeque, front_padding
#import conductor.snbparams
from conductor.vegparams import load_veg_parms

GLACIER_ID = Band.glacier_id
OPEN_GROUND_ID = Band.open_ground_id

cell_ids = ['12345', '23456']

# initially we have just 4 bands loaded for cell '12345', and 3 for cell '23456'
band_ids = {'12345': [0, 1, 2, 3],
            '23456': [1, 2, 3]}
band_size = 100
num_snow_bands = 5

# initial median band elevations
test_median_elevs_simple = [2035, 2172, 2241, 2315]
test_median_elevs = {'12345': [2035, 2172, 2241, 2315],
                    '23456': [1921, 2045, 2164]} 
# initial band map of lower elevation bounds for all existing (valid) bands
test_band_map = {'12345': [2000, 2100, 2200, 2300, 0], # allows for glacier growth at top
                '23456': [0, 1900, 2000, 2100, 0]} # allows for glacier growth at top, and revelation of lower band at bottom

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

num_hrus = {'12345': [2, 3, 2, 1],
            '23456': [3, 3, 2] }


test_veg_types = [11, 19, 22]

test_root_zone_parms = {'11': [0.10, 0.60, 0.20, 0.25, 1.70, 0.15], # 11
                        '19': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0], # 19
                        '22': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0]} # 22


def test_band_simple():
    my_band = Band(test_median_elevs_simple[0])

    assert my_band.median_elev == test_median_elevs_simple[0]
    assert my_band.area_frac == 0
    assert my_band.area_frac_glacier == 0
    assert my_band.area_frac_non_glacier == 0
    assert my_band.area_frac_open_ground == 0
    assert my_band.hrus == {}
    assert my_band.num_hrus == 0

def test_hru_simple():
    my_hru = HydroResponseUnit(test_area_fracs_simple[0], test_root_zone_parms['22'])    

    assert my_hru.area_frac == test_area_fracs_simple[0]
    assert my_hru.root_zone_parms == test_root_zone_parms['22']

def test_band_typical():
    my_band = Band(test_median_elevs_simple[0])

    # Create and populate three HRUs in this Band...
    for veg_type, area_frac, root_zone in zip(test_veg_types, test_area_fracs_simple, test_root_zone_parms):
        my_band.hrus[veg_type] = HydroResponseUnit(area_frac, root_zone)

    assert my_band.median_elev == test_median_elevs_simple[0]
    assert my_band.area_frac == sum(test_area_fracs_simple[0:3])
    assert my_band.area_frac_glacier == test_area_fracs_simple[2]
    assert my_band.area_frac_non_glacier == sum(test_area_fracs_simple[0:3]) - test_area_fracs_simple[2]
    assert my_band.area_frac_open_ground == test_area_fracs_simple[1]
    assert my_band.num_hrus == 3
    for veg, afrac, rzone in zip(test_veg_types, test_area_fracs_simple, test_root_zone_parms):
        assert my_band.hrus[veg].area_frac == afrac
        assert my_band.hrus[veg].root_zone_parms == rzone

# Load up data from large sample vegetation and snow band parameters files and test a few pieces
# of the cells created
def test_merge_cell_input():
    fname = resource_filename('conductor', 'tests/input/snow_band.txt')
    elevation_cells = load_snb_parms(fname, 15)
    fname = resource_filename('conductor', 'tests/input/veg.txt')
    hru_cells = load_veg_parms(fname)
    cells = merge_cell_input(hru_cells, elevation_cells)
    assert len(cells) == 6
    assert len(cells['369560']) == 11
    expected_zs = [ 2076, 2159, 2264, 2354, 2451, 2550, 2620, 2714, 2802 ]
    zs = [ band.median_elev for band in cells['368470'] ]
    assert zs == expected_zs
    expected_afs = {0.000765462339, 0.000873527611, 0.009125511809, 0.009314626034, 0.004426673711, 0.004558753487, 0.001388838859, 0.000737445417}
    afs = { band.hrus[19].area_frac for band in cells['368470'] if 19 in band.hrus }
    assert afs == expected_afs
    assert cells['368470'][0].num_hrus == 2

# Load up the toy problem domain from snow band and vegetation parameter files
# and thoroughly check that the cells were created correctly
def test_cell_simple():
    fname = resource_filename('conductor', 'tests/input/snb_toy_64px.txt')
    elevation_cells = load_snb_parms(fname, num_snow_bands)
    fname = resource_filename('conductor', 'tests/input/vfp_toy_64px.txt')
    hru_cells = load_veg_parms(fname)
    cells = merge_cell_input(hru_cells, elevation_cells)

    # Test that the correct number of Cells were instantiated
    assert len(cells) == 2
    # Test that the correct number of Bands was instantiated for each cell
    assert len(cells[cell_ids[0]]) == 4
    assert len(cells[cell_ids[1]]) == 3
    # Test that the left and right padding is accounted for
    assert cells[cell_ids[0]].left_padding == 0
    assert cells[cell_ids[0]].right_padding == 1
    assert cells[cell_ids[1]].left_padding == 1
    assert cells[cell_ids[1]].right_padding == 1

    # Test that area fractions and root zone parameters for each HRU in each band of one cell are correct
    assert cells['12345'][0].hrus[11].area_frac == test_area_fracs['12345'][0]
    assert cells['12345'][0].hrus[19].area_frac == test_area_fracs['12345'][1]
    assert cells['12345'][1].hrus[11].area_frac == test_area_fracs['12345'][2]
    assert cells['12345'][1].hrus[19].area_frac == test_area_fracs['12345'][3]
    assert cells['12345'][1].hrus[22].area_frac == test_area_fracs['12345'][4]
    assert cells['12345'][2].hrus[19].area_frac == test_area_fracs['12345'][5]
    assert cells['12345'][2].hrus[22].area_frac == test_area_fracs['12345'][6]
    assert cells['12345'][3].hrus[19].area_frac == test_area_fracs['12345'][7]

    assert cells['12345'][0].hrus[11].root_zone_parms == test_root_zone_parms['11']
    assert cells['12345'][0].hrus[19].root_zone_parms == test_root_zone_parms['19']
    assert cells['12345'][1].hrus[11].root_zone_parms == test_root_zone_parms['11']
    assert cells['12345'][1].hrus[19].root_zone_parms == test_root_zone_parms['19']
    assert cells['12345'][1].hrus[22].root_zone_parms == test_root_zone_parms['22']
    assert cells['12345'][2].hrus[19].root_zone_parms == test_root_zone_parms['19']
    assert cells['12345'][2].hrus[22].root_zone_parms == test_root_zone_parms['22']
    assert cells['12345'][3].hrus[19].root_zone_parms == test_root_zone_parms['19']

    for band_id in band_ids['12345']:
        # Test that the number of HRUs reported for each Band in one cell is correct
        assert cells['12345'][band_id].num_hrus == num_hrus['12345'][band_id]
        # Test that all HRU area fractions within a Band add up to original input
        assert sum(hru.area_frac for hru in cells['12345'][band_id].hrus.values()) == sum(test_area_fracs_by_band['12345'][str(band_id)])

def test_cells_dynamic():
    # First, ensure that all Bands' HRU area fractions provided for the test sum to 1
    assert sum(test_area_fracs[cell_ids[0]]) == 1
    assert sum(test_area_fracs[cell_ids[1]]) == 1

    fname = resource_filename('conductor', 'tests/input/snb_toy_64px.txt')
    elevation_cells = load_snb_parms(fname, 5)
    fname = resource_filename('conductor', 'tests/input/vfp_toy_64px.txt')
    hru_cells = load_veg_parms(fname)
    cells = merge_cell_input(hru_cells, elevation_cells)

    ### 2. Simulate glacier expansion over all open ground in Band 2
    new_glacier_area_frac = 0.1875 # 12/64 pixels in toy problem domain, 12/12 pixels for Band 2
    # Glacier HRU area fraction change:
    cells[cell_ids[0]][2].hrus[22].area_frac = new_glacier_area_frac
    # open ground HRU is now gone:
    new_open_ground_area_frac = 0 # not used
    cells[cell_ids[0]][2].delete_hru(OPEN_GROUND_ID)
    # Check that there is only one HRU left in this band
    assert cells[cell_ids[0]][2].num_hrus == 1  

    ### 3. Simulate glacier expansion to replace all open ground in Band 3 (HRU delete + create) 
    new_glacier_area_frac = 0.0625 # 4/64 pixels in toy problem domain. 4/4 in Band 3  
    # open ground HRU is now gone:
    new_open_ground_area_frac = 0
    cells[cell_ids[0]][3].delete_hru(OPEN_GROUND_ID)
    # Confirm that there are (temporarily) no HRUs in this band
    assert cells[cell_ids[0]][3].num_hrus == 0
    # create new glacier HRU:
    cells[cell_ids[0]][3].create_hru(GLACIER_ID, new_glacier_area_frac, test_root_zone_parms['22'])
    # Check that there is only the one glacier HRU in this band
    assert cells[cell_ids[0]][3].num_hrus == 1
    assert cells[cell_ids[0]][3].hrus[22].area_frac == new_glacier_area_frac
    assert cells[cell_ids[0]][3].hrus[22].root_zone_parms == test_root_zone_parms['22']

    ### 4. Simulate glacier growth to create a new elevation Band 4 (stealing one pixel from Band 3)
    new_glacier_area_frac = 0.015625 # 1/64 pixels in domain. 1/1 in Band 4
    # For consistency over whole domain, adjust Band 3 to compensate (this is normally taken care of by update of band_areas):
    cells[cell_ids[0]][3].hrus[22].area_frac -= 0.015625 # area_frac should now be 0.0625 - 0.015625 = 0.046875
    # Confirm existing number of Bands is 4
    assert len(cells[cell_ids[0]]) == 4
    # New Band's initial (single toy pixel) median elevation:
    pixel_elev = 2450
    # Create new Band
    Cell.create_band(cells[cell_ids[0]], pixel_elev)

    # Check that number of Bands has grown by one, and has no HRUs (yet)
    assert len(cells[cell_ids[0]]) == 5
    assert cells[cell_ids[0]][4].num_hrus == 0
    # Create the corresponding new glacier HRU
    cells[cell_ids[0]][4].create_hru(GLACIER_ID, new_glacier_area_frac, test_root_zone_parms['22'])
    # Confirm that this new HRU was correctly instantiated
    assert cells[cell_ids[0]][4].num_hrus == 1
    assert cells[cell_ids[0]][4].hrus[22].area_frac == new_glacier_area_frac
    assert cells[cell_ids[0]][4].hrus[22].root_zone_parms == test_root_zone_parms['22']
    # Confirm this Band's total area_frac is equal to that of its one HRU, and related quantities
    assert cells[cell_ids[0]][4].area_frac == new_glacier_area_frac
    assert cells[cell_ids[0]][4].area_frac_glacier == new_glacier_area_frac
    assert cells[cell_ids[0]][4].area_frac_non_glacier == 0
    assert cells[cell_ids[0]][4].area_frac_open_ground == 0

    ## 5. Simulate an attempt to grow the glacier into a new elevation Band 5 (no 0 pad available)
    pixel_elev = 2550
    with pytest.raises(IndexError):
        Cell.create_band(cells[cell_ids[0]], pixel_elev)
    # Confirm the number of bands has not changed
    assert len(cells[cell_ids[0]]) == 5

    ## 6. Simulate glacier recession completely out of elevation Band 4 (i.e. delete the Band)
    new_glacier_area_frac = 0
    Cell.delete_band(cells[cell_ids[0]], 4)
    # Confirm that there are 4 Bands in total for this cell
    assert len(cells[cell_ids[0]]) == 4
    # For consistency over whole domain, adjust Band 3 to compensate (this is normally taken care of by update of band_areas):
    cells[cell_ids[0]][3].hrus[22].area_frac += 0.015625
    # Confirm that all Band area fractions for this cell still sum to 1
    assert sum(cells[cell_ids[0]][band_idx].area_frac for band_idx, band in enumerate(cells[cell_ids[0]])) == 1

    ## 7. Simulate glacier recession from the lowest existing band, to reveal a yet 
    # lower elevation band (consisting of a single pixel).  This is done in the second test cell, '23456':
    # Initial conditions setup, and check to confirm that no Band 0 currently exists, nor can anything be written to it:
    assert cells[cell_ids[1]].left_padding == 1
    assert cells[cell_ids[1]][0] == None
    with pytest.raises(Exception):
        cells[cell_ids[1]][0].median_elev = 9999
    # New band 0:
    pixel_elev = 1855
    Cell.create_band(cells[cell_ids[1]], pixel_elev)
    # Confirm that the new band was correctly placed in the first position for this cell
    assert cells[cell_ids[1]].left_padding == 0
    assert cells[cell_ids[1]].right_padding == 1
    # Confirm that there are now 4 valid Bands for this cell
    assert len(cells[cell_ids[1]]) == 4
    # Create an open ground HRU in this new lowest band
    new_open_ground_area_frac = 777 # NOTE: this is not a realistic number; just for testing
    cells[cell_ids[1]][0].create_hru(OPEN_GROUND_ID, new_open_ground_area_frac, test_root_zone_parms['19'])
    assert cells[cell_ids[1]][0].num_hrus == 1
    assert cells[cell_ids[1]][0].area_frac == 777

    ## 8. Simulate the glacier re-covering (deleting) the lowest band
    Cell.delete_band(cells[cell_ids[1]], 0)
    # Confirm update of padding
    assert cells[cell_ids[1]].left_padding == 1
    assert cells[cell_ids[1]].right_padding == 1
    # Confirm that there are now 3 valid Bands for this cell
    assert len(cells[cell_ids[1]]) == 3

    ## 9. Attempt to delete a band from the middle of the valid bands
    with pytest.raises(ValueError):
        Cell.delete_band(cells[cell_ids[1]], 2)

    #print('{}'.format(cells))
    #assert 1 == 0 # just to reveal the contents of cells via the print statement above for debugging


def test_update_area_fracs():
    pass

