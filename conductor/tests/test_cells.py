''' This is a set of tests for the cells.py module.
    See conftest.py for details on the test fixtures used.
'''

from collections import OrderedDict
from pkg_resources import resource_filename

import pytest

from conductor.cells import *

GLACIER_ID = Band.glacier_id
OPEN_GROUND_ID = Band.open_ground_id

# Initially we have 4 valid bands loaded for cell 0, and 3 for cell 1
expected_band_ids = {'12345': [0, 1, 2, 3, 4],
                    '23456': [0, 1, 2, 3, 4]}

expected_num_hrus = {'12345': [2, 3, 2, 1, 0],
            '23456': [0, 3, 3, 2, 0] }

expected_root_zone_parms = {'11': [0.10, 0.60, 0.20, 0.25, 1.70, 0.15], # 11
                    '19': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0], # 19
                    '22': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0]} # 22


@pytest.mark.incremental
class TestsSimpleUnit:
    def test_band_and_hru_units(self, simple_unit_test_parms, large_merge_cells_unit_test_parms):

        test_median_elevs_simple, test_median_elevs, test_area_fracs_simple, test_area_fracs, \
            test_area_fracs_by_band, test_veg_types \
            = simple_unit_test_parms

        elevation_cells, hru_cells, expected_zs, expected_afs = large_merge_cells_unit_test_parms

        def test_band_simple(self):
            my_band = Band(test_median_elevs_simple[0])

            assert my_band.median_elev == test_median_elevs_simple[0]
            assert my_band.area_frac == 0
            assert my_band.area_frac_glacier == 0
            assert my_band.area_frac_non_glacier == 0
            assert my_band.area_frac_open_ground == 0
            assert my_band.hrus == {}
            assert my_band.num_hrus == 0

        def test_hru_simple(self):
            my_hru = HydroResponseUnit(test_area_fracs_simple[0], expected_root_zone_parms['22'])

            assert my_hru.area_frac == test_area_fracs_simple[0]
            assert my_hru.root_zone_parms == expected_root_zone_parms['22']

        def test_band_typical(self):
            my_band = Band(test_median_elevs_simple[0])

            # Create and populate three HRUs in this Band...
            for veg_type, area_frac, root_zone in zip(test_veg_types, test_area_fracs_simple, expected_root_zone_parms):
                my_band.hrus[veg_type] = HydroResponseUnit(area_frac, root_zone)

            assert my_band.median_elev == test_median_elevs_simple[0]
            assert my_band.area_frac == sum(test_area_fracs_simple[0:3])
            assert my_band.area_frac_glacier == test_area_fracs_simple[2]
            assert my_band.area_frac_non_glacier == sum(test_area_fracs_simple[0:3]) - test_area_fracs_simple[2]
            assert my_band.area_frac_open_ground == test_area_fracs_simple[1]
            assert my_band.num_hrus == 3
            for veg, afrac, rzone in zip(test_veg_types, test_area_fracs_simple, expected_root_zone_parms):
                assert my_band.hrus[veg].area_frac == afrac
                assert my_band.hrus[veg].root_zone_parms == rzone

       # Load up data from large sample vegetation and snow band parameters files and test a few pieces
        # of the cells created
        def test_merge_cell_input(self):
            cells = merge_cell_input(hru_cells, elevation_cells)
            assert len(cells) == 6
            assert len(cells['369560']) == 15
            zs = [ band.median_elev for band in cells['368470'] ]
            assert zs == expected_zs
            afs = { band.hrus[19].area_frac for band in cells['368470'] if 19 in band.hrus }
            assert afs == expected_afs
            assert cells['368470'][0].num_hrus == 2

        test_band_simple(self)
        test_hru_simple(self)
        test_band_typical(self)
        test_merge_cell_input(self)

    def test_cell_creation_simple(self, simple_unit_test_parms, toy_domain_64px_cells):
        """Load up the toy problem domain from snow band and vegetation parameter files
        and thoroughly check (against simple unit test parms, which should also partially
        represent the toy domain) that the cells were created correctly"""

        test_median_elevs_simple, test_median_elevs, test_area_fracs_simple, test_area_fracs, \
            test_area_fracs_by_band, test_veg_types = simple_unit_test_parms

        cells, cell_ids, num_snow_bands, band_size, cellid_map, surf_dem, glacier_mask, \
            cell_band_pixel_elevations = toy_domain_64px_cells

        # Test that the correct number of Cells were instantiated
        assert len(cells) == len(cell_ids)
        # Test that the correct number of Bands was instantiated for each cell (including dummy bands)
        assert len(cells[cell_ids[0]]) == num_snow_bands
        assert len(cells[cell_ids[1]]) == num_snow_bands
        # Test that the number of HRUs in the valid bands are correct (including dummy bands)
        assert cells[cell_ids[0]][0].num_hrus == expected_num_hrus[cell_ids[0]][0]
        assert cells[cell_ids[0]][1].num_hrus == expected_num_hrus[cell_ids[0]][1]
        assert cells[cell_ids[0]][2].num_hrus == expected_num_hrus[cell_ids[0]][2]
        assert cells[cell_ids[0]][3].num_hrus == expected_num_hrus[cell_ids[0]][3]
        assert cells[cell_ids[0]][4].num_hrus == expected_num_hrus[cell_ids[0]][4]

        assert cells[cell_ids[1]][0].num_hrus == expected_num_hrus[cell_ids[1]][0]
        assert cells[cell_ids[1]][1].num_hrus == expected_num_hrus[cell_ids[1]][1]
        assert cells[cell_ids[1]][2].num_hrus == expected_num_hrus[cell_ids[1]][2]
        assert cells[cell_ids[1]][3].num_hrus == expected_num_hrus[cell_ids[1]][3]
        assert cells[cell_ids[1]][4].num_hrus == expected_num_hrus[cell_ids[1]][4]

        # Test that area fractions and root zone parameters for each HRU in each band of one cell are correct
        assert cells[cell_ids[0]][0].hrus[11].area_frac == test_area_fracs[cell_ids[0]][0]
        assert cells[cell_ids[0]][0].hrus[19].area_frac == test_area_fracs[cell_ids[0]][1]
        assert cells[cell_ids[0]][1].hrus[11].area_frac == test_area_fracs[cell_ids[0]][2]
        assert cells[cell_ids[0]][1].hrus[19].area_frac == test_area_fracs[cell_ids[0]][3]
        assert cells[cell_ids[0]][1].hrus[22].area_frac == test_area_fracs[cell_ids[0]][4]
        assert cells[cell_ids[0]][2].hrus[19].area_frac == test_area_fracs[cell_ids[0]][5]
        assert cells[cell_ids[0]][2].hrus[22].area_frac == test_area_fracs[cell_ids[0]][6]
        assert cells[cell_ids[0]][3].hrus[19].area_frac == test_area_fracs[cell_ids[0]][7]

        assert cells[cell_ids[0]][0].hrus[11].root_zone_parms == expected_root_zone_parms['11']
        assert cells[cell_ids[0]][0].hrus[19].root_zone_parms == expected_root_zone_parms['19']
        assert cells[cell_ids[0]][1].hrus[11].root_zone_parms == expected_root_zone_parms['11']
        assert cells[cell_ids[0]][1].hrus[19].root_zone_parms == expected_root_zone_parms['19']
        assert cells[cell_ids[0]][1].hrus[22].root_zone_parms == expected_root_zone_parms['22']
        assert cells[cell_ids[0]][2].hrus[19].root_zone_parms == expected_root_zone_parms['19']
        assert cells[cell_ids[0]][2].hrus[22].root_zone_parms == expected_root_zone_parms['22']
        assert cells[cell_ids[0]][3].hrus[19].root_zone_parms == expected_root_zone_parms['19']

        for band_id in expected_band_ids[cell_ids[0]]:
            # Test that the number of HRUs reported for each Band in one cell is correct
            assert cells[cell_ids[0]][band_id].num_hrus == expected_num_hrus[cell_ids[0]][band_id]
            # Test that all HRU area fractions within a Band add up to original input
            assert sum(hru.area_frac for hru in cells[cell_ids[0]][band_id].hrus.values()) == sum(test_area_fracs_by_band[cell_ids[0]][str(band_id)])

@pytest.mark.incremental
class TestsDynamic:
    def test_cells_dynamic(self, toy_domain_64px_cells):

        cells, cell_ids, num_snow_bands, band_size, cellid_map, surf_dem, glacier_mask, \
            cell_band_pixel_elevations = toy_domain_64px_cells

        def test_existing_glacier_growth_within_band_replacing_all_open_ground(self):
            """test_cells_dynamic -- Test #1: Simulates glacier expansion over all open ground in Band 2 """
            new_glacier_area_frac = 0.1875 # 12/64 pixels in toy problem domain, 12/12 pixels for Band 2
            # Glacier HRU area fraction change:
            cells[cell_ids[0]][2].hrus[22].area_frac = new_glacier_area_frac
            # open ground HRU is now gone:
            new_open_ground_area_frac = 0 # not used
            cells[cell_ids[0]][2].delete_hru(OPEN_GROUND_ID)
            # Check that there is only one HRU left in this band
            assert cells[cell_ids[0]][2].num_hrus == 1

        def test_new_glacier_growth_into_band_and_replacing_all_open_ground(self):
            """test_cells_dynamic -- Test #2: Simulates glacier expansion to replace all open ground in Band 3 (HRU delete + create)"""
            new_glacier_area_frac = 0.0625 # 4/64 pixels in toy problem domain. 4/4 in Band 3  
            # open ground HRU is now gone:
            new_open_ground_area_frac = 0
            cells[cell_ids[0]][3].delete_hru(OPEN_GROUND_ID)
            # Confirm that there are (temporarily) no HRUs in this band
            assert cells[cell_ids[0]][3].num_hrus == 0
            # create new glacier HRU:
            cells[cell_ids[0]][3].create_hru(GLACIER_ID, new_glacier_area_frac)
            # Check that there is only the one glacier HRU in this band
            assert cells[cell_ids[0]][3].num_hrus == 1
            assert cells[cell_ids[0]][3].hrus[22].area_frac == new_glacier_area_frac
            assert cells[cell_ids[0]][3].hrus[22].root_zone_parms == expected_root_zone_parms['22']


        def test_new_glacier_growth_into_upper_dummy_band(self):
            """test_cells_dynamic -- Test #3: Simulates glacier growth to create a new elevation Band 4 (stealing one pixel from Band 3)"""
            new_glacier_area_frac = 0.015625 # 1/64 pixels in domain. 1/1 in Band 4
            # For consistency over whole domain, adjust Band 3 to compensate (this is normally taken care of by update of band_areas):
            cells[cell_ids[0]][3].hrus[22].area_frac -= 0.015625 # area_frac should now be 0.0625 - 0.015625 = 0.046875
            # Confirm that there are currently no HRUs in the new band         
            assert cells[cell_ids[0]][4].num_hrus == 0
            # Create the corresponding new glacier HRU
            cells[cell_ids[0]][4].create_hru(GLACIER_ID, new_glacier_area_frac)
            # Confirm that this new HRU was correctly instantiated
            assert cells[cell_ids[0]][4].num_hrus == 1
            assert cells[cell_ids[0]][4].hrus[22].area_frac == new_glacier_area_frac
            assert cells[cell_ids[0]][4].hrus[22].root_zone_parms == expected_root_zone_parms['22']
            # Confirm this Band's total area_frac is equal to that of its one HRU, and related quantities
            assert cells[cell_ids[0]][4].area_frac == new_glacier_area_frac
            assert cells[cell_ids[0]][4].area_frac_glacier == new_glacier_area_frac
            assert cells[cell_ids[0]][4].area_frac_non_glacier == 0
            assert cells[cell_ids[0]][4].area_frac_open_ground == 0

        test_existing_glacier_growth_within_band_replacing_all_open_ground(self)
        test_new_glacier_growth_into_band_and_replacing_all_open_ground(self)
        test_new_glacier_growth_into_upper_dummy_band(self)

@pytest.mark.incremental
class TestsAreaFracUpdate:
    def test_update_area_fracs(self, toy_domain_64px_cells, toy_domain_64px_rgm_vic_map_file_readout):
        cells, cell_ids, num_snow_bands, band_size, cellid_map, initial_surf_dem, initial_glacier_mask, \
            cell_band_pixel_elevations = toy_domain_64px_cells

        _, _, cell_areas, num_cols_dem, num_rows_dem = toy_domain_64px_rgm_vic_map_file_readout

        def test_no_changes(self):
            cells_orig = deepcopy(cells)

            update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
                      initial_surf_dem, num_rows_dem, num_cols_dem, initial_glacier_mask)

            assert cells == cells_orig

        test_no_changes(self)

        def test_glacier_growth_over_open_ground_and_vegetation_in_band(self):
            """ Simulates Band 1 losing some of its open ground area to glacier growth
                (confined to cell '12345' here) """
            surf_dem = deepcopy(initial_surf_dem)
            # TODO: modify surf_dem accordingly, and run update_area_fracs()


        def test_glacier_growth_over_open_ground_and_vegetation_in_band(self):
            """ Simulates Band 1 losing all its open ground and some vegetated area to glacier growth"""

        def test_glacier_growth_over_remaining_vegetation_in_band(self):
            """ Simulates Band 1 losing (some of) its only remaining non-glacier (vegetated) HRU to glacier growth"""
            pass

        def test_attempt_new_glacier_growth_into_unavailable_lower_band(self):
            """test_cells_dynamic -- Test #7: Simulates the glacier expanding downward into an elevation band
            that has not been anticipated, i.e. not enough 0 pads were provided on the lower end in the snow band file"""
            pass

        def test_glacier_growth_to_conceal_lowest_band(self):
            """test_cells_dynamic -- Test #8: Simulates the glacier re-covering the lowest band thickly enough
            such that the pixels elevations in that area no longer belong to that band (i.e. the band must be deleted).
            NOTE: the glacier area fraction for existing Band 1 is not being updated in this test"""
            pass

        def test_existing_glacier_shrink_revealing_new_lower_band(self):
            """test_cells_dynamic -- Test #6: Simulates glacier recession from the lowest existing band, to reveal a yet
            lower elevation band (consisting of a single pixel).  This is done in the second test cell, ID '23456' """
            pass

        def test_attempt_new_glacier_growth_into_unavailable_higher_band(self):
            """test_cells_dynamic -- Test #4: Simulates a (failing) attempt to grow the glacier into a new elevation Band 5 (no 0 pad available)"""
            pass
