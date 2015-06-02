''' This is a set of tests for the cells.py module
'''

import pytest

from collections import OrderedDict

#import conductor.cells
from conductor.cells import *

GLACIER_ID = 22
OPEN_GROUND_ID = 19

cell_ids = ['12345', '23456']
band_ids = ['0', '1', '2', '3']
band_size = 100

test_median_elevs = [2050, 2150, 2250, 2350]
test_band_map = {'12345': [2000, 2100, 2200, 2300, 0], '23456': [1900, 2000, 2100, 0, 0]}

test_area_fracs = [0.1875, 0.25, # Band 0 (11, 19)
                    0.0625, 0.125, 0.125, # Band 1 (11, 19, 22)
                    0.0625, 0.125, # Band 2 (19, 22)
                    0.0625] # Band 3 (19) 
test_veg_types = [11, 19, 22]
test_root_zone_parms = [[0.10, 0.60, 0.20, 0.25, 1.70, 0.15], # 11
                        [0.1, 1.0, 0.1, 0.0, 0.1, 0.0], # 19
                        [0.1, 1.0, 0.1, 0.0, 0.1, 0.0]] # 22

def test_band_simple():
    my_band = Band(test_median_elevs[0], GLACIER_ID, OPEN_GROUND_ID)

    assert my_band.median_elev == test_median_elevs[0]
    assert my_band.area_frac == 0
    assert my_band.area_frac_glacier == 0
    assert my_band.area_frac_non_glacier == 0
    assert my_band.area_frac_open_ground == 0
    assert my_band.hrus == []
    assert my_band.num_hrus == 0

def test_hru_simple():
    my_hru = HydroResponseUnit(GLACIER_ID, test_area_fracs[0], test_root_zone_parms[2])    

    assert my_hru.veg_type == GLACIER_ID
    assert my_hru.area_frac == test_area_fracs[0]
    assert my_hru.root_zone_parms == test_root_zone_parms[2]

def test_band_typical():
    my_band = Band(test_median_elevs[0], GLACIER_ID, OPEN_GROUND_ID)

    # Create and populate three HRUs in this Band...
    # Tree HRU:
    my_band.hrus.append(HydroResponseUnit(test_veg_types[0], test_area_fracs[0], test_root_zone_parms[0]))
    # Open ground HRU:
    my_band.hrus.append(HydroResponseUnit(test_veg_types[1], test_area_fracs[1], test_root_zone_parms[1]))
    # Glacier HRU:
    my_band.hrus.append(HydroResponseUnit(test_veg_types[2], test_area_fracs[2], test_root_zone_parms[2]))

    assert my_band.median_elev == test_median_elevs[0]
    assert my_band.area_frac == sum(test_area_fracs[0:3])
    assert my_band.area_frac_glacier == test_area_fracs[2]
    assert my_band.area_frac_non_glacier == sum(test_area_fracs[0:3]) - test_area_fracs[2]
    assert my_band.area_frac_open_ground == test_area_fracs[1]
    assert my_band.num_hrus == 3
    assert my_band.hrus[0].veg_type == test_veg_types[0]
    assert my_band.hrus[0].area_frac == test_area_fracs[0]
    assert my_band.hrus[0].root_zone_parms == test_root_zone_parms[0]
    assert my_band.hrus[1].veg_type == test_veg_types[1]
    assert my_band.hrus[1].area_frac == test_area_fracs[1]
    assert my_band.hrus[1].root_zone_parms == test_root_zone_parms[1]
    assert my_band.hrus[2].veg_type == test_veg_types[2]
    assert my_band.hrus[2].area_frac == test_area_fracs[2]
    assert my_band.hrus[2].root_zone_parms == test_root_zone_parms[2]

def test_cells_simple():
    cells = OrderedDict()

    # Create a cell with 3 Bands
    cells[cell_ids[0]] = OrderedDict()
    cells[cell_ids[0]][band_ids[0]] = Band(test_median_elevs[0], GLACIER_ID, OPEN_GROUND_ID)
    # Create and populate three HRUs in this Band...
    # Tree HRU:
    cells[cell_ids[0]][band_ids[0]].hrus.append(HydroResponseUnit(test_veg_types[0], test_area_fracs[0], test_root_zone_parms[0]))
    # Open ground HRU:
    cells[cell_ids[0]][band_ids[0]].hrus.append(HydroResponseUnit(test_veg_types[1], test_area_fracs[1], test_root_zone_parms[1]))
    # Glacier HRU:
    cells[cell_ids[0]][band_ids[0]].hrus.append(HydroResponseUnit(test_veg_types[2], test_area_fracs[2], test_root_zone_parms[2]))

    cells[cell_ids[0]][band_ids[1]] = Band(test_median_elevs[1], GLACIER_ID, OPEN_GROUND_ID)
    cells[cell_ids[0]][band_ids[2]] = Band(test_median_elevs[2], GLACIER_ID, OPEN_GROUND_ID)
    cells[cell_ids[0]][band_ids[3]] = Band(test_median_elevs[3], GLACIER_ID, OPEN_GROUND_ID)

    # Create another cell with 3 Bands
    cells[cell_ids[1]] = OrderedDict()
    cells[cell_ids[1]][band_ids[0]] = Band(test_median_elevs[0], GLACIER_ID, OPEN_GROUND_ID)
    cells[cell_ids[1]][band_ids[1]] = Band(test_median_elevs[1], GLACIER_ID, OPEN_GROUND_ID)
    cells[cell_ids[1]][band_ids[2]] = Band(test_median_elevs[2], GLACIER_ID, OPEN_GROUND_ID)
    cells[cell_ids[1]][band_ids[3]] = Band(test_median_elevs[3], GLACIER_ID, OPEN_GROUND_ID)

    # Test that all HRU area fractions within a Band add up to original input
    assert sum(hru.area_frac for hru in cells[cell_ids[0]][band_ids[0]].hrus) == sum(test_area_fracs[0:3])

def test_cells_dynamic():

    # Ensure that all Bands' HRU area fractions provided for test sum to 1
    assert sum(test_area_fracs) == 1

    cells = OrderedDict()

    ## 1. Set up initial conditions
    # Create a cell with 3 Bands
    cells[cell_ids[0]] = OrderedDict()
    # Band 0:
    cells[cell_ids[0]][band_ids[0]] = Band(test_median_elevs[0], GLACIER_ID, OPEN_GROUND_ID)
    # Tree HRU:
    cells[cell_ids[0]][band_ids[0]].hrus.append(HydroResponseUnit(test_veg_types[0], test_area_fracs[0], test_root_zone_parms[0]))
    # Open ground HRU:
    cells[cell_ids[0]][band_ids[0]].hrus.append(HydroResponseUnit(test_veg_types[1], test_area_fracs[1], test_root_zone_parms[1]))

    # Band 1:
    cells[cell_ids[0]][band_ids[1]] = Band(test_median_elevs[1], GLACIER_ID, OPEN_GROUND_ID)
    # Tree HRU:
    cells[cell_ids[0]][band_ids[1]].hrus.append(HydroResponseUnit(test_veg_types[0], test_area_fracs[2], test_root_zone_parms[0]))
    # Open ground HRU:
    cells[cell_ids[0]][band_ids[1]].hrus.append(HydroResponseUnit(test_veg_types[1], test_area_fracs[3], test_root_zone_parms[1]))
    # Glacier HRU:
    cells[cell_ids[0]][band_ids[1]].hrus.append(HydroResponseUnit(test_veg_types[2], test_area_fracs[4], test_root_zone_parms[2]))

    # Band 2:
    cells[cell_ids[0]][band_ids[2]] = Band(test_median_elevs[2], GLACIER_ID, OPEN_GROUND_ID)
    # Open ground HRU:
    cells[cell_ids[0]][band_ids[2]].hrus.append(HydroResponseUnit(test_veg_types[1], test_area_fracs[5], test_root_zone_parms[1]))
    # Glacier HRU:
    cells[cell_ids[0]][band_ids[2]].hrus.append(HydroResponseUnit(test_veg_types[2], test_area_fracs[6], test_root_zone_parms[2]))

    # Band 3:
    cells[cell_ids[0]][band_ids[3]] = Band(test_median_elevs[3], GLACIER_ID, OPEN_GROUND_ID)
    # Open ground HRU:
    cells[cell_ids[0]][band_ids[3]].hrus.append(HydroResponseUnit(test_veg_types[1], test_area_fracs[7], test_root_zone_parms[1]))

    ### 2. Simulate glacier expansion over all open ground in Band 2
    new_glacier_area_frac = 0.1875 # 12/64 pixels in toy problem domain, 12/12 pixels for Band 2
    # Glacier HRU area fraction change:
    cells[cell_ids[0]][band_ids[2]].hrus[1].area_frac = new_glacier_area_frac
    # open ground HRU is now gone:
    new_open_ground_area_frac = 0
    delete_hru(cells, cell_ids[0], band_ids[2], OPEN_GROUND_ID)
    # Check that there is only one HRU left in this band
    assert cells[cell_ids[0]][band_ids[3]].num_hrus == 1    

    ### 3. Simulate glacier expansion to replace all open ground in Band 3 (HRU delete + create) 
    new_glacier_area_frac = 0.0625 # 4/64 pixels in toy problem domain. 4/4 in Band 3  
    # open ground HRU is now gone:
    new_open_ground_area_frac = 0
    delete_hru(cells, cell_ids[0], band_ids[3], OPEN_GROUND_ID)
    # create new glacier HRU:
    create_hru(cells, cell_ids[0], band_ids[3], GLACIER_ID, new_glacier_area_frac, test_root_zone_parms[2])
    # assign area fraction to this new glacier HRU:
    cells[cell_ids[0]][band_ids[3]].hrus[0].area_frac = new_glacier_area_frac
    # Check that there is only one HRU left in this band
    assert cells[cell_ids[0]][band_ids[3]].num_hrus == 1
    assert cells[cell_ids[0]][band_ids[3]].hrus[0].veg_type == GLACIER_ID
    assert cells[cell_ids[0]][band_ids[3]].hrus[0].area_frac == new_glacier_area_frac
    assert cells[cell_ids[0]][band_ids[3]].hrus[0].root_zone_parms == test_root_zone_parms[2]

    ### 4. Simulate glacier growth to create a new elevation Band 4 (stealing one pixel from Band 3)
    new_glacier_area_frac = 0.015625 # 1/64 pixels in domain. 1/1 in Band 4
    # For consistency over whole domain, adjust Band 3 to compensate (this is normally taken care of by update of band_areas):
    cells[cell_ids[0]][band_ids[3]].hrus[0].area_frac -= 0.015625
    # New Band's initial (single toy pixel) median elevation:
    pixel_elev = 2450
    new_band_idx = create_band(cells, cell_ids[0], pixel_elev, band_size, test_band_map, GLACIER_ID, OPEN_GROUND_ID)
    new_band_id = str(new_band_idx) # this is true because we only ever append to upper end of band_map
#    print('new_band_id: {}'.format(new_band_id))
#    print('band_map: {}'.format(test_band_map))
    # Check that this new band was created in the band_map at new_band_idx
    assert test_band_map[cell_ids[0]][new_band_idx] == 2400
    # Test that band_map stayed the same length (num_snow_bands = 5 in this case)
    assert len(test_band_map[cell_ids[0]]) == 5
    # Create the corresponding new glacier HRU
    create_hru(cells, cell_ids[0], new_band_id, GLACIER_ID, new_glacier_area_frac, test_root_zone_parms[2])
    # Check out this new HRU
    assert cells[cell_ids[0]][new_band_id].num_hrus == 1
    assert cells[cell_ids[0]][new_band_id].hrus[0].veg_type == GLACIER_ID
    assert cells[cell_ids[0]][new_band_id].hrus[0].area_frac == new_glacier_area_frac
    assert cells[cell_ids[0]][new_band_id].hrus[0].root_zone_parms == test_root_zone_parms[2]
    # Confirm this Band's total area_frac is equal to that of its one HRU, and related quantities
    assert cells[cell_ids[0]][new_band_id].area_frac == new_glacier_area_frac
    assert cells[cell_ids[0]][new_band_id].area_frac_glacier == new_glacier_area_frac
    assert cells[cell_ids[0]][new_band_id].area_frac_non_glacier == 0
    assert cells[cell_ids[0]][new_band_id].area_frac_open_ground == 0

    ## Simulate an attempt to grow the glacier into a new elevation Band 5 (no 0 pad available)
    pixel_elev = 2550
    with pytest.raises(NameError):
        new_band_idx = create_band(cells, cell_ids[0], pixel_elev, band_size, test_band_map, GLACIER_ID, OPEN_GROUND_ID)

    ## 5. Simulate glacier recession completely out of elevation Band 4 (i.e. delete the Band)
    new_glacier_area_frac = 0
    delete_band(cells, cell_ids[0], 2400, test_band_map)
    # Confirm that there are 4 Bands in total for this cell
    assert len(cells[cell_ids[0]]) == 4
    # Confirm that the entry for this former valid band is set to 0 in band_map
    assert test_band_map[cell_ids[0]][4] == 0
    # For consistency over whole domain, adjust Band 3 to compensate (this is normally taken care of by update of band_areas):
    cells[cell_ids[0]][band_ids[3]].hrus[0].area_frac += 0.015625

    ## 6. Confirm that all Band area fractions for this cell still sum to 1
    assert sum(cells[cell_ids[0]][band].area_frac for band in cells[cell_ids[0]]) == 1