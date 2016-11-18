from pkg_resources import resource_filename

from conductor.vegparams import load_veg_parms

def test_load_veg_parms():
  fname = resource_filename('conductor', 'tests/input/veg.txt')
  cells = load_veg_parms(fname)
  assert len(cells) == 6
  assert len(cells['368470']) == 16
