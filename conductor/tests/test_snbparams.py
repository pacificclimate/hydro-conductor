from pkg_resources import resource_filename

import pytest

from conductor.snbparams import load_snb_parms

def test_load_snb_parms():
  fname = resource_filename('conductor', 'tests/input/snow_band.txt')
  cells = load_snb_parms(fname, 15)
  assert len(cells) == 6
  assert len(cells['369560']) == 15
  expected_zs = [ 2076, 2159, 2264, 2354, 2451, 2550, 2620, 2714, 2802,\
    2950, 3050, 3150, 3250, 3350, 3450 ]
  zs = [ band.median_elev for band in cells['368470'] ]
  assert zs == expected_zs

