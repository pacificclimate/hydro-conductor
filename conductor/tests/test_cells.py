''' This is a set of tests for the cells.py module.
  See conftest.py for details on the test fixtures used.
'''

from collections import OrderedDict
from pkg_resources import resource_filename
from copy import deepcopy

import pytest
import mock

from conductor.cells import *
from conductor.io import update_glacier_mask

GLACIER_ID = Band.glacier_id
OPEN_GROUND_ID = Band.open_ground_id

# Initially we have 4 valid bands loaded for cell 0, and 3 for cell 1
expected_band_ids = {
  '12345': [0, 1, 2, 3, 4],
  '23456': [0, 1, 2, 3, 4]
}

expected_num_hrus = {
  '12345': [2, 3, 2, 1, 0],
  '23456': [0, 3, 3, 2, 0]
}

expected_root_zone_parms = {
  '11': [0.10, 0.60, 0.20, 0.25, 1.70, 0.15], # 11
  '19': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0], # 19
  '22': [0.1, 1.0, 0.1, 0.0, 0.1, 0.0] # 22
}

@pytest.mark.incremental
class TestsSimpleUnit:
  def test_band_and_hru_units(self, simple_unit_test_parms,\
    large_merge_cells_unit_test_parms):

    test_median_elevs_simple, test_median_elevs, test_area_fracs_simple,\
      test_area_fracs, test_area_fracs_by_band, test_veg_types \
      = simple_unit_test_parms

    elevation_cells, hru_cells, expected_zs, expected_afs \
      = large_merge_cells_unit_test_parms

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
      my_hru = HydroResponseUnit(test_area_fracs_simple[0],\
        expected_root_zone_parms['22'])

      assert my_hru.area_frac == test_area_fracs_simple[0]
      assert my_hru.root_zone_parms == expected_root_zone_parms['22']

    def test_band_typical(self):
      my_band = Band(test_median_elevs_simple[0])

      # Create and populate three HRUs in this Band...
      for veg_type, area_frac, root_zone in zip(test_veg_types,\
        test_area_fracs_simple, expected_root_zone_parms):
          my_band.hrus[veg_type] = HydroResponseUnit(area_frac, root_zone)

      assert my_band.median_elev == test_median_elevs_simple[0]
      assert my_band.area_frac == sum(test_area_fracs_simple[0:3])
      assert my_band.area_frac_glacier == test_area_fracs_simple[2]
      assert my_band.area_frac_non_glacier == sum(test_area_fracs_simple[0:3])\
        - test_area_fracs_simple[2]
      assert my_band.area_frac_open_ground == test_area_fracs_simple[1]
      assert my_band.num_hrus == 3
      for veg, afrac, rzone in zip(test_veg_types, test_area_fracs_simple,\
        expected_root_zone_parms):
        assert my_band.hrus[veg].area_frac == afrac
        assert my_band.hrus[veg].root_zone_parms == rzone

    # Load up data from large sample vegetation and snow band parameters files
    # and test a few pieces of the cells created
    def test_merge_cell_input(self):
      cells = merge_cell_input(hru_cells, elevation_cells)
      assert len(cells) == 6
      assert len(cells['369560'].bands) == 15
      zs = [ band.median_elev for band in cells['368470'].bands ]
      assert zs == expected_zs
      afs = { band.hrus[19].area_frac for band in cells['368470'].bands \
        if 19 in band.hrus }
      assert afs == expected_afs
      assert cells['368470'].bands[0].num_hrus == 2

    test_band_simple(self)
    test_hru_simple(self)
    test_band_typical(self)
    test_merge_cell_input(self)

  def test_cell_creation_simple(self, simple_unit_test_parms,\
    toy_domain_64px_cells):
    """Load up the toy problem domain from snow band and vegetation parameter
    files and thoroughly check (against simple unit test parms, which should
    also partially represent the toy domain) that the cells were created
    correctly
    """

    test_median_elevs_simple, test_median_elevs, test_area_fracs_simple,\
      test_area_fracs, test_area_fracs_by_band, test_veg_types\
      = simple_unit_test_parms

    cells, cell_ids, num_snow_bands, band_size, cellid_map, bed_dem,\
      surf_dem, glacier_mask, cell_band_pixel_elevations\
      = toy_domain_64px_cells

    # Test that the correct number of Cells were instantiated
    assert len(cells) == len(cell_ids)
    # Test that the correct number of Bands was instantiated for each cell 
    # (including dummy bands)
    assert len(cells[cell_ids[0]].bands) == num_snow_bands
    assert len(cells[cell_ids[1]].bands) == num_snow_bands
    # Test that the number of HRUs in the valid bands are correct (including
    # dummy bands)
    assert cells[cell_ids[0]].bands[0].num_hrus == expected_num_hrus[cell_ids[0]][0]
    assert cells[cell_ids[0]].bands[1].num_hrus == expected_num_hrus[cell_ids[0]][1]
    assert cells[cell_ids[0]].bands[2].num_hrus == expected_num_hrus[cell_ids[0]][2]
    assert cells[cell_ids[0]].bands[3].num_hrus == expected_num_hrus[cell_ids[0]][3]
    assert cells[cell_ids[0]].bands[4].num_hrus == expected_num_hrus[cell_ids[0]][4]

    assert cells[cell_ids[1]].bands[0].num_hrus == expected_num_hrus[cell_ids[1]][0]
    assert cells[cell_ids[1]].bands[1].num_hrus == expected_num_hrus[cell_ids[1]][1]
    assert cells[cell_ids[1]].bands[2].num_hrus == expected_num_hrus[cell_ids[1]][2]
    assert cells[cell_ids[1]].bands[3].num_hrus == expected_num_hrus[cell_ids[1]][3]
    assert cells[cell_ids[1]].bands[4].num_hrus == expected_num_hrus[cell_ids[1]][4]

    # Test that area fractions and root zone parameters for each HRU in each 
    # band of one cell are correct
    assert cells[cell_ids[0]].bands[0].hrus[11].area_frac == test_area_fracs[cell_ids[0]][0]
    assert cells[cell_ids[0]].bands[0].hrus[19].area_frac == test_area_fracs[cell_ids[0]][1]
    assert cells[cell_ids[0]].bands[1].hrus[11].area_frac == test_area_fracs[cell_ids[0]][2]
    assert cells[cell_ids[0]].bands[1].hrus[19].area_frac == test_area_fracs[cell_ids[0]][3]
    assert cells[cell_ids[0]].bands[1].hrus[22].area_frac == test_area_fracs[cell_ids[0]][4]
    assert cells[cell_ids[0]].bands[2].hrus[19].area_frac == test_area_fracs[cell_ids[0]][5]
    assert cells[cell_ids[0]].bands[2].hrus[22].area_frac == test_area_fracs[cell_ids[0]][6]
    assert cells[cell_ids[0]].bands[3].hrus[19].area_frac == test_area_fracs[cell_ids[0]][7]

    assert cells[cell_ids[0]].bands[0].hrus[11].root_zone_parms == expected_root_zone_parms['11']
    assert cells[cell_ids[0]].bands[0].hrus[19].root_zone_parms == expected_root_zone_parms['19']
    assert cells[cell_ids[0]].bands[1].hrus[11].root_zone_parms == expected_root_zone_parms['11']
    assert cells[cell_ids[0]].bands[1].hrus[19].root_zone_parms == expected_root_zone_parms['19']
    assert cells[cell_ids[0]].bands[1].hrus[22].root_zone_parms == expected_root_zone_parms['22']
    assert cells[cell_ids[0]].bands[2].hrus[19].root_zone_parms == expected_root_zone_parms['19']
    assert cells[cell_ids[0]].bands[2].hrus[22].root_zone_parms == expected_root_zone_parms['22']
    assert cells[cell_ids[0]].bands[3].hrus[19].root_zone_parms == expected_root_zone_parms['19']

    for band_id in expected_band_ids[cell_ids[0]]:
      # Test that the number of HRUs reported for each Band in a cell is correct
      assert cells[cell_ids[0]].bands[band_id].num_hrus\
        == expected_num_hrus[cell_ids[0]][band_id]
      # Test that all HRU area fractions within a Band add up to original input
      assert sum(hru.area_frac for hru in cells[cell_ids[0]].bands[band_id].hrus.values())\
        == sum(test_area_fracs_by_band[cell_ids[0]][str(band_id)])

@pytest.mark.incremental
class TestsDynamic:
  def test_cells_dynamic(self, toy_domain_64px_cells):

    cells, cell_ids, num_snow_bands, band_size, cellid_map, bed_dem, surf_dem,\
      glacier_mask, cell_band_pixel_elevations = toy_domain_64px_cells

    def test_existing_glacier_growth_within_band_replacing_all_open_ground(self):
      """test_cells_dynamic -- Test #1: Simulates glacier expansion over all
        open ground in Band 2
      """
      new_glacier_area_frac = 0.1875 # 12/64 pixels in toy problem domain, 12/12 pixels for Band 2
      # Glacier HRU area fraction change:
      cells[cell_ids[0]].bands[2].hrus[22].area_frac = new_glacier_area_frac
      # open ground HRU is now gone:
      new_open_ground_area_frac = 0 # not used
      cells[cell_ids[0]].bands[2].delete_hru(OPEN_GROUND_ID)
      # Check that there is only one HRU left in this band
      assert cells[cell_ids[0]].bands[2].num_hrus == 1

    def test_new_glacier_growth_into_band_and_replacing_all_open_ground(self):
      """test_cells_dynamic -- Test #2: Simulates glacier expansion to replace
        all open ground in Band 3 (HRU delete + create)
        """
      # 4/64 pixels in toy problem domain. 4/4 in Band 3:
      new_glacier_area_frac = 0.0625 
      # open ground HRU is now gone:
      new_open_ground_area_frac = 0
      cells[cell_ids[0]].bands[3].delete_hru(OPEN_GROUND_ID)
      # Confirm that there are (temporarily) no HRUs in this band
      assert cells[cell_ids[0]].bands[3].num_hrus == 0
      # create new glacier HRU:
      cells[cell_ids[0]].bands[3].create_hru(GLACIER_ID, new_glacier_area_frac)
      # Check that there is only the one glacier HRU in this band
      assert cells[cell_ids[0]].bands[3].num_hrus == 1
      assert cells[cell_ids[0]].bands[3].hrus[22].area_frac == new_glacier_area_frac
      assert cells[cell_ids[0]].bands[3].hrus[22].root_zone_parms\
        == expected_root_zone_parms['22']

    def test_new_glacier_growth_into_upper_dummy_band(self):
      """test_cells_dynamic -- Test #3: Simulates glacier growth to create a
        new elevation Band 4 (stealing one pixel from Band 3)
      """
      new_glacier_area_frac = 0.015625 # 1/64 pixels in domain. 1/1 in Band 4
      # For consistency over whole domain, adjust Band 3 to compensate
      # (this is normally taken care of by update of band_areas):
      cells[cell_ids[0]].bands[3].hrus[22].area_frac -= 0.015625
      # (area_frac should now be 0.0625 - 0.015625 = 0.046875)

      # Confirm that there are currently no HRUs in the new band         
      assert cells[cell_ids[0]].bands[4].num_hrus == 0
      # Create the corresponding new glacier HRU
      cells[cell_ids[0]].bands[4].create_hru(GLACIER_ID, new_glacier_area_frac)
      # Confirm that this new HRU was correctly instantiated
      assert cells[cell_ids[0]].bands[4].num_hrus == 1
      assert cells[cell_ids[0]].bands[4].hrus[22].area_frac == new_glacier_area_frac
      assert cells[cell_ids[0]].bands[4].hrus[22].root_zone_parms\
        == expected_root_zone_parms['22']
      # Confirm this Band's total area_frac is equal to that of its one HRU, 
      # and related quantities
      assert cells[cell_ids[0]].bands[4].area_frac == new_glacier_area_frac
      assert cells[cell_ids[0]].bands[4].area_frac_glacier == new_glacier_area_frac
      assert cells[cell_ids[0]].bands[4].area_frac_non_glacier == 0
      assert cells[cell_ids[0]].bands[4].area_frac_open_ground == 0

    test_existing_glacier_growth_within_band_replacing_all_open_ground(self)
    test_new_glacier_growth_into_band_and_replacing_all_open_ground(self)
    test_new_glacier_growth_into_upper_dummy_band(self)

def mock_update_hru_state(source_hru, dest_hru, case, new_area_fracs):
  """ Mock function for cells.update_hru_state(), just returns the
    state update case that was given in the case input parameter.
  """
  return case

@pytest.mark.incremental
class TestsAreaFracUpdate:
  def test_update_area_fracs(self, toy_domain_64px_cells,\
    toy_domain_64px_rgm_vic_map_file_readout):
    cells, cell_ids, num_snow_bands, band_size, cellid_map, bed_dem,\
    initial_surf_dem, initial_glacier_mask, cell_band_pixel_elevations\
      = toy_domain_64px_cells

    _, _, cell_areas, num_cols_dem, num_rows_dem\
      = toy_domain_64px_rgm_vic_map_file_readout

    # We'll use this copy to make modifications to as the DEM evolves through
    # the incremental tests
    surf_dem = deepcopy(initial_surf_dem)
    # thickness of np.NaN pads around DEM, glacier_mask, and cellid_map
    dem_padding_thickness = 2

    def test_no_changes(self):
      cells_orig = deepcopy(cells)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        initial_surf_dem, num_rows_dem, num_cols_dem, initial_glacier_mask)

      assert cells == cells_orig

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_growth_over_some_open_ground_in_band(self, mock_update_hru_state_fcn):
      """ Simulates Band 2 of cell '12345' losing some of its open ground
        area to glacier growth.
        
        Should trigger state update CASE 3 (glacier HRU expansion)
        and CASE 3 again (open ground HRU shrinkage).

        Initial surf_dem for cell '12345':
        [
          [2065, 2055, 2045, 2035, 2025, 2015, 2005, 2000],
          [2075, 2100, 2120, 2140, 2130, 2120, 2100, 2005],
          [2085, 2110, 2250, 2270, 2260, 2240, 2110, 2010],
          [2090, 2120, 2260, 2377, 2310, 2250, 2125, 2015],
          [2070, 2110, 2250, 2340, 2320, 2250, 2130, 2020],
          [2090, 2105, 2200, 2210, 2220, 2220, 2120, 2015],
          [2090, 2100, 2105, 2110, 2140, 2150, 2130, 2010],
          [2080, 2075, 2065, 2055, 2045, 2035, 2020, 2000]
        ]

        Changed elevations due to glacier growth:
        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, 2230, 2240, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
        """
      surf_dem[dem_padding_thickness + 5][dem_padding_thickness + 2 :\
        dem_padding_thickness + 4] = [2230, 2240]

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3 twice
      assert mock_update_hru_state_fcn.call_count == 2
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[1]
      assert args[2] == '3'

      assert cells['12345'].bands[2].num_hrus == 2
      assert cells['12345'].bands[2].area_frac == 12/64
      assert cells['12345'].bands[2].area_frac_open_ground == 2/64
      assert cells['12345'].bands[2].area_frac_glacier == 10/64

      # Total number of valid bands
      assert len([band for band in cells['12345'].bands if band.num_hrus > 0]) == 4

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_growth_over_remaining_open_ground_in_band(self, mock_update_hru_state_fcn):
      """ Simulates Band 2 of cell '12345' losing all its remaining open ground. 

        This should trigger state update CASE 3 (glacier expansion) and
        CASE 4b (loss of remaining open ground in band).

        Changed elevations due to glacier growth (incremental from last test):
        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, 2240, 2230, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 5][dem_padding_thickness + 4 :\
        dem_padding_thickness + 6] = [2240, 2230]

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3 and CASE 4b
      assert mock_update_hru_state_fcn.call_count == 2
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[1]
      assert args[2] == '4b'

      assert cells['12345'].bands[2].num_hrus == 1
      assert cells['12345'].bands[2].area_frac == 12/64
      assert cells['12345'].bands[2].area_frac_open_ground == 0
      assert cells['12345'].bands[2].area_frac_glacier == 12/64

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_growth_over_remaining_open_ground_and_some_vegetation_in_band(self, mock_update_hru_state_fcn):
      """ Simulates Band 1 of cell '12345' losing all its open ground and some 
        vegetated area to glacier growth.

        This should trigger state update CASE 3 (glacier expansion) and
        CASE 4b (loss of remaining open ground in band) and CASE 3 again
        (loss of some vegetation).

        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, 2120, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, 2130, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, 2145, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, 2150, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, 2140, xxxx],
          [xxxx, xxxx, xxxx, 2140, 2160, 2160, 2150, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 1][dem_padding_thickness + 6] = 2120
      surf_dem[dem_padding_thickness + 2][dem_padding_thickness + 6] = 2130
      surf_dem[dem_padding_thickness + 3][dem_padding_thickness + 6] = 2145
      surf_dem[dem_padding_thickness + 4][dem_padding_thickness + 6] = 2150
      surf_dem[dem_padding_thickness + 5][dem_padding_thickness + 6] = 2140
      surf_dem[dem_padding_thickness + 6][dem_padding_thickness + 3:\
        dem_padding_thickness + 7] = [2140, 2160, 2160, 2150]

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3, 4b, and 3
      assert mock_update_hru_state_fcn.call_count == 3
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[1]
      assert args[2] == '4b'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[2]
      assert args[2] == '3'

      assert cells['12345'].bands[1].num_hrus == 2
      assert cells['12345'].bands[1].area_frac == 20/64
      assert cells['12345'].bands[1].area_frac_open_ground == 0
      assert cells['12345'].bands[1].area_frac_glacier == 17/64
      assert cells['12345'].bands[1].hrus[11].area_frac == 3/64

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_growth_over_remaining_vegetation_in_band(self, mock_update_hru_state_fcn):
      """ Simulates Band 1 of cell '12345' losing its remaining vegetated
        HRU to glacier growth.

        This should trigger state update CASE 3 (glacier expansion) and
        CASE 4b (loss of remaining vegetation in band).

        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, 2115, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, 2110, 2125, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 5][dem_padding_thickness + 1] = 2115
      surf_dem[dem_padding_thickness + 6][dem_padding_thickness + 1 :\
        dem_padding_thickness + 3] = [2110, 2125]

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3 and 4b
      assert mock_update_hru_state_fcn.call_count == 2
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[1]
      assert args[2] == '4b'

      assert cells['12345'].bands[1].num_hrus == 1
      assert cells['12345'].bands[1].area_frac == 20/64
      assert cells['12345'].bands[1].area_frac_open_ground == 0
      assert cells['12345'].bands[1].area_frac_glacier == 20/64

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_growth_into_band_with_no_existing_glacier(self, mock_update_hru_state_fcn):
      """ Simulates Band 0 of cell '12345' acquiring a new glacier HRU.

        This should trigger state update CASE 1 (trivial - a new open
        ground HRU is created but update_hru_state is not called) and
        CASE 3 (loss of some open ground).

        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, 2030],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, 2040],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 7] = 2030
      surf_dem[dem_padding_thickness + 1][dem_padding_thickness + 7] = 2040

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3
      assert mock_update_hru_state_fcn.call_count == 1
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'

      assert cells['12345'].bands[0].num_hrus == 3
      assert cells['12345'].bands[0].area_frac == 28/64
      assert cells['12345'].bands[0].area_frac_open_ground == 14/64
      assert cells['12345'].bands[0].area_frac_glacier == 2/64
      assert cells['12345'].bands[0].hrus[11].area_frac == 12/64

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_receding_to_reveal_open_ground_in_band(self, mock_update_hru_state_fcn):
      """ Simulates Band 1 of cell '12345', which is completely covered in
        glacier, ceding some area to open ground.

        This should trigger state update CASE 3 (glacier expansion)
        and CASE 1 (trivial - a new open ground HRU is created but
        update_hru_state is not called).

        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, 2100, 2105, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 6][dem_padding_thickness + 1 :\
        dem_padding_thickness + 3] = [2100, 2105]

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3
      assert mock_update_hru_state_fcn.call_count == 1
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'

      assert cells['12345'].bands[1].num_hrus == 2
      assert cells['12345'].bands[1].area_frac == 20/64
      assert cells['12345'].bands[1].area_frac_open_ground == 2/64
      assert cells['12345'].bands[1].area_frac_glacier == 18/64

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_receding_further_in_band(self, mock_update_hru_state_fcn):
      """ Simulates Band 1 of cell '12345' ceding additional area to open
        ground (2 pixels which were open and tree types at the very beginning)

        This should trigger state update CASE 3 (glacier shrink) and
        CASE 3 again (open ground expansion).

        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, 2110, 2140, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 6][dem_padding_thickness + 3 :\
        dem_padding_thickness + 5] = [2110, 2140]

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3 twice
      assert mock_update_hru_state_fcn.call_count == 2
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[1]
      assert args[2] == '3'

      assert cells['12345'].bands[1].num_hrus == 2
      assert cells['12345'].bands[1].area_frac == 20/64
      assert cells['12345'].bands[1].area_frac_open_ground == 4/64
      assert cells['12345'].bands[1].area_frac_glacier == 16/64

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_existing_glacier_shrink_revealing_new_lower_band(self, mock_update_hru_state_fcn):
      """ Simulates glacier recession out of the lowest existing band of cell
        '23456', to reveal a yet lower elevation band (consisting of one pixel).

        Available bands as per snow band parameter file:
        '12345': [2000, 2100, 2200, 2300, 0] # allows for glacier growth at top
        '23456': [0, 1900, 2000, 2100, 0] # allows for glacier growth at top,
        and revelation of lower band at bottom

        This should trigger state update CASE 3 (glacier shrink) and
        CASE 1 (trivial - open ground HRU creation at newly revealed band).

        [
          [xxxx, xxxx, 1850, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      # Total number of valid bands before
      assert len([band for band in cells['23456'].bands if band.num_hrus > 0]) == 3

      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 8 + 2] = 1850

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3
      assert mock_update_hru_state_fcn.call_count == 1
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'

      assert cells['23456'].bands[1].num_hrus == 3
      assert cells['23456'].bands[1].area_frac == 27/64
      assert cells['23456'].bands[1].area_frac_open_ground == 10/64
      assert cells['23456'].bands[1].area_frac_glacier == 1/64
      assert cells['23456'].bands[1].hrus[11].area_frac == 16/64

      # New lowest band
      assert cells['23456'].bands[0].num_hrus == 1
      assert cells['23456'].bands[0].area_frac == 1/64
      assert cells['23456'].bands[0].area_frac_open_ground == 1/64
      assert cells['23456'].bands[0].area_frac_glacier == 0

      # Total number of valid bands after
      assert len([band for band in cells['23456'].bands if band.num_hrus > 0]) == 4

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_growth_into_new_lower_band(self, mock_update_hru_state_fcn):
      """ Simulates glacier growing back over the pixel of the new lowest band
        in cell '23456' (from the previous test), but at a lesser thickness
        such that the pixel is still within Band 0.

        This should trigger state update CASE 1 (trivial - glacier creation
        at the lowest band) and CASE 4b (loss of all open ground HRU at the
        lowest band).

        [
          [xxxx, xxxx, 1880, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      # Copy initial lowest band's area_frac, to use for re-initializing at end
      initial_area_frac = cells['23456'].bands[0].area_frac
      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 8 + 2] = 1880

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 4b
      assert mock_update_hru_state_fcn.call_count == 1
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '4b'

      assert cells['23456'].bands[0].num_hrus == 1
      assert cells['23456'].bands[0].area_frac == 1/64
      assert cells['23456'].bands[0].area_frac_open_ground == 0
      assert cells['23456'].bands[0].area_frac_glacier == 1/64

      # Reinstate original elevation of changed pixel and remove created
      # glacier HRU for next incremental test 
      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 8 + 2] = 1850
      cells['23456'].bands[0].delete_hru(22)
      cells['23456'].bands[0].create_hru(19, initial_area_frac)

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_thickening_to_conceal_lowest_band_of_open_ground(self, mock_update_hru_state_fcn):
      """ Simulates the glacier growing over open ground areas lying in the
        new lowest band of cell '23456' so thick that the pixels elevations
        in that area no longer belong to that band (i.e. all HRUs in the band
        must be deleted).

        This should trigger state update CASE 3 (glacier in band 1 expanding)
        and CASE 5d (where both the glacier and band area fractions in Band 0
        become zero, and there is no lower band to transfer state to).

        [
          [xxxx, xxxx, 1905, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 8 + 2] = 1905

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3 and 5d
      assert mock_update_hru_state_fcn.call_count == 2
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[1]
      assert args[2] == '5d'

      assert cells['23456'].bands[0].num_hrus == 0 # we delete open ground HRUs
      assert cells['23456'].bands[0].lower_bound == 1800
      assert cells['23456'].bands[0].median_elev == 1800
      assert cells['23456'].bands[0].area_frac == 0
      assert cells['23456'].bands[0].area_frac_open_ground == 0
      assert cells['23456'].bands[0].area_frac_glacier == 0

      assert cells['23456'].bands[1].num_hrus == 3
      assert cells['23456'].bands[1].area_frac == 28/64
      assert cells['23456'].bands[1].area_frac_open_ground == 10/64
      assert cells['23456'].bands[1].area_frac_glacier == 2/64
      assert cells['23456'].bands[1].hrus[11].area_frac == 16/64

      # Total number of valid bands after (should not include the lowest one
      # now, because HRU was deleted)
      assert len([band for band in cells['23456'].bands if band.num_hrus > 0]) == 3

      # Reinstate lowest band with single glaciated pixel at 1880m for the next test
      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 8 + 2] = 1880
      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)
      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_growth_into_new_higher_band(self, mock_update_hru_state_fcn):
      """ Simulates glacier growing in thickness from the highest existing
        valid band in cell '23456' into a new higher Band 4 (for which there
        is a 0 pad in the snow band file to accommodate it). 

        This should trigger state update CASE 3 (glacier in Band 3 shrinking)
        and CASE 1 (trivial - no call to update_hru_state) where a new glacier
        HRU is created in Band 4.

        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, 2200, 2210, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 3][dem_padding_thickness + 8 + 3 :\
        dem_padding_thickness + 8 + 5] = [2200, 2210]

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3
      assert mock_update_hru_state_fcn.call_count == 1
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'

      assert cells['23456'].bands[3].num_hrus == 2
      assert cells['23456'].bands[3].area_frac == 14/64
      assert cells['23456'].bands[3].area_frac_open_ground == 8/64
      assert cells['23456'].bands[3].area_frac_glacier == 6/64

      assert cells['23456'].bands[4].num_hrus == 1
      assert cells['23456'].bands[4].area_frac == 2/64
      assert cells['23456'].bands[4].area_frac_open_ground == 0
      assert cells['23456'].bands[4].area_frac_glacier == 2/64

      # Total number of valid bands after
      assert len([band for band in cells['23456'].bands if band.num_hrus > 0]) == 5

    def test_attempt_new_glacier_shrink_into_unavailable_lower_band(self):
      """ Simulates a (failing) attempt to grow the glacier into a new yet
        lower elevation band (where there is no 0 pad available in the snow
        band parameter file) in cell '23456'
        [
          [xxxx, xxxx, xxxx, 1799, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 8 + 3] = 1799

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      with pytest.raises(Exception) as message:
        update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
          surf_dem, num_rows_dem, num_cols_dem, glacier_mask)
      assert 'One or more RGM output DEM pixels lies below the bounds of the \
lowest defined elevation band (< 1800.0m) as defined by the Snow Band \
Parameter File for cell 23456. You may need to add or shift the zero \
padding to accommodate this.' in str(message.value)

      # Remove error condition by reinstating glacier in offending pixel in
      # the surface DEM for the next test
      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 8 + 3] = 1995
      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)
      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

    def test_attempt_new_glacier_growth_into_unavailable_higher_band(self):
      """ Simulates a (failing) attempt to grow the glacier into a new yet
        higher elevation band (where there is no 0 pad available in the snow
        band parameter file) in cell '23456' 
        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, 2300, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 3][dem_padding_thickness + 8 + 3] = 2300

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      with pytest.raises(Exception) as message:
        update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
          surf_dem, num_rows_dem, num_cols_dem, glacier_mask)
      assert 'One or more RGM output DEM pixels lies above the bounds of \
the highest defined elevation band (>= 2300.0m) as defined by the \
Snow Band Parameter File for cell 23456. You may need to add or shift \
the zero padding to accommodate this.' in str(message.value)

       # Remove error condition by reinstating glacier in offending pixel in
       # the surface DEM for the next test
      surf_dem[dem_padding_thickness + 3][dem_padding_thickness + 8 + 3] = 2200
      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)
      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_thickening_to_conceal_lowest_band_of_glacier(self, mock_update_hru_state_fcn):
      """ Simulates the glacier thickening over areas lying in the lowest
        band of cell '23456' so much that the pixels elevations in that area
        no longer belong to that band (i.e. all HRUs in the band must be
        deleted, except glacier which is set to zero area fraction).

        This should trigger state update CASE 3 (glacier in Band 1 expanding)
        and CASE 5d (where both the glacier and band area fractions in Band 0
        become zero, and there is no lower band to transfer state to).

        [
          [xxxx, xxxx, 1900, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 0][dem_padding_thickness + 8 + 2] = 1900

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem,\
        num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
        surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 3 and 5d
      assert mock_update_hru_state_fcn.call_count == 2
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '3'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[1]
      assert args[2] == '5d'

      # we never delete glacier HRUs, vis-a-vis VIC's shadow glaciers:
      assert cells['23456'].bands[0].num_hrus == 1
      assert cells['23456'].bands[0].lower_bound == 1800
      assert cells['23456'].bands[0].median_elev == 1800
      assert cells['23456'].bands[0].area_frac == 0
      assert cells['23456'].bands[0].area_frac_open_ground == 0
      assert cells['23456'].bands[0].area_frac_glacier == 0
 
      assert cells['23456'].bands[1].num_hrus == 3
      assert cells['23456'].bands[1].area_frac == 28/64
      assert cells['23456'].bands[1].area_frac_open_ground == 10/64
      assert cells['23456'].bands[1].area_frac_glacier == 2/64
      assert cells['23456'].bands[1].hrus[11].area_frac == 16/64

      # Total number of valid bands after (should include the lowest one now,
      # because of the glacier HRU)
      assert len([band for band in cells['23456'].bands if band.num_hrus > 0]) == 5

    @mock.patch('conductor.cells.update_hru_state', side_effect=mock_update_hru_state)
    def test_glacier_receding_from_top_band_leaving_band_area_as_zero(self, mock_update_hru_state_fcn):
      """ Simulates the glacier receding out of the highest band of cell
        '23456' entirely, which consisted only of glacier HRUs, thus
        leaving that band's area fraction as zero.

        This should trigger state update CASE 5a (glacier in Band 4
        disappears and state is transferred to glacier in Band 3)
        and CASE 3 (glacier in Band 3 expands).

        [
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, 2150, 2170, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx],
          [xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx, xxxx]
        ]
      """
      surf_dem[dem_padding_thickness + 3][dem_padding_thickness + 8 + 3 : dem_padding_thickness + 8 + 5] = [2150, 2170]

      glacier_mask = update_glacier_mask(surf_dem, bed_dem, num_rows_dem, num_cols_dem)

      update_area_fracs(cells, cell_areas, cellid_map, num_snow_bands,\
                surf_dem, num_rows_dem, num_cols_dem, glacier_mask)

      # Confirm that update_hru_state() was called for CASE 5a and 3
      assert mock_update_hru_state_fcn.call_count == 2
      args, kwargs = mock_update_hru_state_fcn.call_args_list[0]
      assert args[2] == '5a'
      args, kwargs = mock_update_hru_state_fcn.call_args_list[1]
      assert args[2] == '3'

      assert cells['23456'].bands[4].num_hrus == 1 # shadow glacier HRU remains
      assert cells['23456'].bands[4].area_frac == 0
      assert cells['23456'].bands[4].area_frac_open_ground == 0
      assert cells['23456'].bands[4].area_frac_glacier == 0

      assert cells['23456'].bands[3].num_hrus == 2
      assert cells['23456'].bands[3].area_frac == 16/64
      assert cells['23456'].bands[3].area_frac_open_ground == 8/64
      assert cells['23456'].bands[3].area_frac_glacier == 8/64

      # Total number of valid bands after
      assert len([band for band in cells['23456'].bands if band.num_hrus > 0]) == 5

    def test_confirm_final_state(self):
      """ Final test to confirm that the final state of both grid cells is as
        expected after the sequence of operations performed upon them in the
        preceding tests.
      """
      assert cells['12345'].bands[0].num_hrus == 3
      assert cells['12345'].bands[0].lower_bound == 2000
      assert cells['12345'].bands[0].median_elev == \
        np.median([2065, 2055, 2045, 2035, 2025, 2015, 2005, 2030, 2040, 2010,\
          2015, 2020, 2015, 2010, 2000, 2080, 2075, 2065, 2055, 2045, 2035,\
          2020, 2075, 2085, 2090, 2070, 2090, 2090])
      assert cells['12345'].bands[0].area_frac == 28/64
      assert cells['12345'].bands[0].area_frac_open_ground == 14/64
      assert cells['12345'].bands[0].area_frac_glacier == 2/64
      assert cells['12345'].bands[0].hrus[11].area_frac == 12/64

      assert cells['12345'].bands[1].num_hrus == 2
      assert cells['12345'].bands[1].lower_bound == 2100
      assert cells['12345'].bands[1].median_elev == \
        np.median([2100, 2120, 2140, 2130, 2120, 2120, 2130, 2145, 2150, 2140,\
          2150, 2100, 2105, 2110, 2140, 2160, 2110, 2120, 2120, 2115])
      assert cells['12345'].bands[1].area_frac == 20/64
      assert cells['12345'].bands[1].area_frac_open_ground == 4/64
      assert cells['12345'].bands[1].area_frac_glacier == 16/64

      assert cells['12345'].bands[2].num_hrus == 1
      assert cells['12345'].bands[2].lower_bound == 2200
      assert cells['12345'].bands[2].median_elev == \
        np.median([2250, 2270, 2260, 2240, 2250, 2250, 2230, 2230, 2240, 2240,\
          2260, 2250])
      assert cells['12345'].bands[2].area_frac == 12/64
      assert cells['12345'].bands[2].area_frac_open_ground == 0/64
      assert cells['12345'].bands[2].area_frac_glacier == 12/64

      assert cells['12345'].bands[3].num_hrus == 1
      assert cells['12345'].bands[3].lower_bound == 2300
      assert cells['12345'].bands[3].median_elev == np.median([2377, 2310, 2340, 2320])
      assert cells['12345'].bands[3].area_frac == 4/64
      assert cells['12345'].bands[3].area_frac_open_ground == 4/64
      assert cells['12345'].bands[3].area_frac_glacier == 0/64

      assert len([band for band in cells['12345'].bands if band.num_hrus > 0]) == 4

      # the shadow glacier HRU for the currently invalid lowest band:
      assert cells['23456'].bands[0].num_hrus == 1
      assert cells['23456'].bands[0].lower_bound == 1800
      assert cells['23456'].bands[0].median_elev == 1800
      assert cells['23456'].bands[0].area_frac == 0/64
      assert cells['23456'].bands[0].area_frac_open_ground == 0/64
      assert cells['23456'].bands[0].area_frac_glacier == 0/64

      assert cells['23456'].bands[1].num_hrus == 3
      assert cells['23456'].bands[1].lower_bound == 1900
      assert cells['23456'].bands[1].median_elev == \
        np.median([1970, 1975, 1900, 1995, 1975, 1965, 1960, 1960, 1965, 1970,\
          1975, 1980, 1980, 1970, 1960, 1965, 1965, 1970, 1970, 1975, 1960,\
          1950, 1970, 1975, 1985, 1990, 1980, 1970])
      assert cells['23456'].bands[1].area_frac == 28/64
      assert cells['23456'].bands[1].area_frac_open_ground == 10/64
      assert cells['23456'].bands[1].area_frac_glacier == 2/64
      assert cells['23456'].bands[1].hrus[11].area_frac == 16/64

      assert cells['23456'].bands[2].num_hrus == 3
      assert cells['23456'].bands[2].lower_bound == 2000
      assert cells['23456'].bands[2].median_elev == \
        np.median([2000, 2045, 2055, 2005, 2005, 2000, 2000, 2000, 2005, 2000,\
          2000, 2000, 2000, 2020, 2035, 2025, 2000, 2005, 2010, 2005])
      assert cells['23456'].bands[2].area_frac == 20/64
      assert cells['23456'].bands[2].area_frac_open_ground == 8/64
      assert cells['23456'].bands[2].area_frac_glacier == 2/64
      assert cells['23456'].bands[2].hrus[11].area_frac == 10/64

      assert cells['23456'].bands[3].num_hrus == 2
      assert cells['23456'].bands[3].lower_bound == 2100
      assert cells['23456'].bands[3].median_elev == np.median([2100, 2155, 2160, 2140,\
                                                         2105, 2150, 2170, 2130,\
                                                         2110, 2150, 2140, 2105,\
                                                         2105, 2105, 2110, 2100])
      assert cells['23456'].bands[3].area_frac == 16/64
      assert cells['23456'].bands[3].area_frac_open_ground == 8/64
      assert cells['23456'].bands[3].area_frac_glacier == 8/64

      assert cells['23456'].bands[4].num_hrus == 1 # shadow glacier HRU
      assert cells['23456'].bands[4].lower_bound == 2200
      assert cells['23456'].bands[4].median_elev == 2200
      assert cells['23456'].bands[4].area_frac == 0/64
      assert cells['23456'].bands[4].area_frac_open_ground == 0/64
      assert cells['23456'].bands[4].area_frac_glacier == 0/64

      # This should include all bands, including the lowest and highest,
      # now each consisting of a shadow glacier HRU
      assert len([band for band in cells['23456'].bands if band.num_hrus > 0]) == 5

    test_no_changes(self)
    test_glacier_growth_over_some_open_ground_in_band(self)
    test_glacier_growth_over_remaining_open_ground_in_band(self)
    test_glacier_growth_over_remaining_open_ground_and_some_vegetation_in_band(self)
    test_glacier_growth_over_remaining_vegetation_in_band(self)
    test_glacier_growth_into_band_with_no_existing_glacier(self)
    test_glacier_receding_to_reveal_open_ground_in_band(self)
    test_glacier_receding_further_in_band(self)

    test_existing_glacier_shrink_revealing_new_lower_band(self)
    test_glacier_growth_into_new_lower_band(self)
    test_glacier_thickening_to_conceal_lowest_band_of_open_ground(self)
    test_glacier_growth_into_new_higher_band(self)
    test_attempt_new_glacier_shrink_into_unavailable_lower_band(self)
    test_attempt_new_glacier_growth_into_unavailable_higher_band(self)
    test_glacier_thickening_to_conceal_lowest_band_of_glacier(self)
    test_glacier_receding_from_top_band_leaving_band_area_as_zero(self)
    test_confirm_final_state(self)
