from pkg_resources import resource_filename

import pytest

from conductor.snbparams import front_padding, load_snb_parms

def test_load_snb_parms():
    fname = resource_filename('conductor', 'tests/input/snow_band.txt')
    cells = load_snb_parms(fname, 15)
    assert len(cells) == 6
    assert len(cells['369560']) == 15
    expected_zs = [ 2076, 2159, 2264, 2354, 2451, 2550, 2620, 2714, 2802, 2900, 3000, 3100, 3200, 3300, 3400 ]
    zs = [ band.median_elev for band in cells['368470'] ]
    assert zs == expected_zs

