import io
from pkg_resources import resource_stream, resource_filename

import pytest

from conductor.cells import *
from conductor.snbparams import load_snb_parms, PaddedDeque, front_padding
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
    stream = resource_stream('conductor', 'tests/input/global.txt')
    return io.TextIOWrapper(stream)

@pytest.fixture(scope="module")
def toy_domain_64px_cells():
    fname = resource_filename('conductor', 'tests/input/snb_toy_64px.txt')
    elevation_cells = load_snb_parms(fname, 5)
    fname = resource_filename('conductor', 'tests/input/vfp_toy_64px.txt')
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
    return cells, cell_ids, num_snow_bands, band_size, expected_band_ids, expected_root_zone_parms

@pytest.fixture(scope="module")
def simple_unit_test_parms():

    # NOT USED BUT USEFUL INFO: initial band map of lower elevation bounds for all existing (valid) bands
    # test_band_map = {'12345': [2000, 2100, 2200, 2300, 0], # allows for glacier growth at top
    #            '23456': [0, 1900, 2000, 2100, 0]} # allows for glacier growth at top, and revelation of lower band at bottom

    # initial median band elevations
    test_median_elevs_simple = [2035, 2172, 2241, 2315]
    test_median_elevs = {'12345': [2035, 2172, 2241, 2315],
                        '23456': [1921, 2045, 2164]} 

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