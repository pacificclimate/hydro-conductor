from pkg_resources import resource_filename

import pytest

from conductor.snbparams import front_padding, load_snb_parms, PaddedDeque

@pytest.mark.parametrize(('input_', 'expected'),
                         [([0, 0, 1, 0], 2),
                          ([1, 1, 1, 1], 0),
                          ([0, 0, 0, 0], 4),
                          ([0, 0, 1, 1], 2),
                          ([1, 1, 0, 0], 0)
                          ])
def test_front_padding(input_, expected):
    assert front_padding(input_) == expected

def test_load_snb_parms():
    fname = resource_filename('conductor', 'tests/input/snow_band.txt')
    cells = load_snb_parms(fname, 15)
    assert len(cells) == 6
    assert len(cells['369560']) == 11
    # the last 4 elevations are dummy bands
    expected_zs = [ 2076, 2159, 2264, 2354, 2451, 2550, 2620, 2714, 2802 ]    
    zs = [ band.median_elev for band in cells['368470'] ]
    assert zs == expected_zs

@pytest.mark.parametrize(('args', 'kwargs', 'padded_result'),
                         # no padding
                         [(([1, 2, 3], 3), {}, [1, 2, 3]),
                          # left padding
                          (([1, 2, 3], 4), {'left_padding': 1}, [None, 1, 2, 3]),
                          # right padding
                          (([1, 2, 3], 4), {}, [1, 2, 3, None]), 
                          (([1, 2, 3], 5), {'left_padding': 1}, [None, 1, 2, 3, None]), # double padded
                          ])
def test_padded_deque_init(args, kwargs, padded_result):
    pd = PaddedDeque(*args, **kwargs)
    assert list(pd.padded_iter()) == padded_result
    assert list(pd) == [1, 2, 3]

def test_padded_deque_empty_insert():
    pd = PaddedDeque([], 3)
    pd.append(1)
    assert list(pd.padded_iter()) == [1, None, None]
    assert list(pd) == [1]

def test_padded_deque_overflow():
    pd = PaddedDeque([1], 1)
    with pytest.raises(IndexError):
        pd.appendleft(4)
    with pytest.raises(IndexError):
        pd.append(4)
    
def test_padded_deque_getitem():
    pd = PaddedDeque([1, 2, 3], 5, left_padding=1)
    assert pd[0] is None
    assert pd[2] == 2
    assert pd[4] is None
    with pytest.raises(IndexError):
        pd[5]

def test_padded_deque_peek():
    pd = PaddedDeque([1, 2, 3], 5, left_padding=1)
    assert pd.peekleft() == 1
    assert pd.peekright() == 3

def test_padded_deque_pop():
    pd = PaddedDeque(range(5), 5)
    assert list(pd) == [0, 1, 2, 3, 4]
    assert pd.popleft() == 0
    assert pd[0] is None
    assert pd.pop() == 4
    assert pd[4] is None
    assert list(pd.padded_iter()) == [None, 1, 2, 3, None]
    assert list(pd) == [1, 2, 3]
    for _ in range(3):
        pd.pop()
    assert list(pd.padded_iter()) == [None] * 5
    assert list(pd) == []
    with pytest.raises(IndexError):
        pd.pop()
