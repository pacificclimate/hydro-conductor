import datetime
import numpy as np

import pytest

from vic_rgm_conductor import run_ranges

@pytest.mark.parametrize(('start', 'end', 'glac', 'expected'), [
  # Simulation ends after the end of water year
  ('1950/01/01', '1959/12/31', '1955/10/01', [
    ('1950/01/01', '1956/09/30'),
    ('1956/10/01', '1957/09/30'),
    ('1957/10/01', '1958/09/30'),
    ('1958/10/01', '1959/09/30'),
  ]),
  # Simulation ends before the end of the water year
  ('1950/01/01', '1959/06/01', '1955/10/01', [
    ('1950/01/01', '1956/09/30'),
    ('1956/10/01', '1957/09/30'),
    ('1957/10/01', '1958/09/30'),
    ('1958/10/01', '1959/06/01'),
  ]),
  # (Almost) No spin up period for VIC
  ('1950/01/01', '1952/12/31', '1950/10/01', [
    ('1950/01/01', '1951/09/30'),
    ('1951/10/01', '1952/09/30'),
  ]),
])
def test_run_ranges(start, end, glac, expected):
  def todate(datestr):
    return datetime.datetime.strptime(datestr, '%Y/%m/%d').date()

  iterator = run_ranges(*[todate(d) for d in [start, end, glac]])
  for (s, e), (exp_s, exp_e) in zip(iterator, expected):
    assert s == todate(exp_s)
    assert e == todate(exp_e)

def test_run_ranges_warning(recwarn):
  """ test a glacier start date not aligned to the water year
  """
  i = run_ranges(None, None, datetime.date(2000, 1, 1))
  next(i) # The warning doesn't appear until we draw from the iterator
  w = recwarn.pop()
  assert 'run_ranges assumes that glacier_start' in str(w.message)

