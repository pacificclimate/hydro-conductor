from conductor.vic_globals import Scalar, Boolean, Filename, Mapping, List, Global

import pytest

def test_scalar_int():
    class Foo(object):
        x = Scalar(int)
    f = Foo()
    f.x = "5"
    assert f.x == 5

    with pytest.raises(ValueError):
        f.x = "10.01"
    
    with pytest.raises(ValueError):
        f.x = "this not convertable to an int"

def test_scalar_float():
    class Bar(object):
        x = Scalar(float)
    b = Bar()
    b.x = "10.125"
    assert b.x == 10.125
    
    with pytest.raises(ValueError):
        b.x = "not a float"

def test_scalar_str():
    class Foo(object):
        x = Scalar(str)
    f = Foo()
    f.x = "10.1"
    assert f.x == "10.1"
        
@pytest.mark.parametrize(('input_', 'expected'), [
    ("TRUE", True),
    ("True", True),
    ("true", True),
    ("FALSE", False),
    ("False", False),
    ("false", False),
    ("0", False),
    ])
def test_boolean(input_, expected):
    class Blah(object):
        x = Boolean()
    b = Blah()
    b.x = input_
    assert b.x == expected

def test_filename():
    class Stuff(object):
        x = Filename()
    s = Stuff()
    s.x = '/tmp/my_file.txt'
    assert s.x == '/tmp/my_file.txt'

    with pytest.raises(ValueError):
        s.x = '/this/directory/does/not/exist'

def test_mapping():
    class Foo(object):
        x = Mapping()
    f = Foo()
    f.x = "some_key some_value"
    assert f.x == {'some_key': 'some_value'}
    f.x = "another_key another value with multiple splits"
    assert f.x == {'some_key': 'some_value',
                   'another_key': 'another value with multiple splits'}

    with pytest.raises(ValueError) as ex:
        f.x = 'Cannot_split_this'
    assert 'requires a key and value separated by whitespace' in str(ex.value)
    #print(Foo.__dict__['x'].__str__(f, Foo, 'x'))

def test_list():
    class Foo(object):
        x = List()
    f = Foo()
    assert f.x == []
    f.x = "Something"
    f.x = 1.0
    f.x = 8
    assert f.x == ["Something", 1.0, 8]

def test_global_init(sample_global_file_string):
    g = Global(sample_global_file_string)
    print(str(g))

def test_force_types(sample_global_file_string):
    g = Global(sample_global_file_string)
    force_type = g._str_member('force_type')
    force_dt = g._str_member('force_dt')

    assert force_dt == 'FORCE_DT 1\n'

    expected_force_types = ['FORCE_TYPE SHORTWAVE SHORTWAVE\n', 'FORCE_TYPE LONGWAVE LONGWAVE\n', \
                        'FORCE_TYPE AIR_TEMP AIR_TEMP\n', 'FORCE_TYPE PRESSURE PRESSURE\n', \
                        'FORCE_TYPE DENSITY DENSITY\n', 'FORCE_TYPE VP VP\n', 'FORCE_TYPE WIND WIND\n', \
                        'FORCE_TYPE PREC PREC\n']
    for ft in expected_force_types:
        assert ft in force_type

# def test_global_parms_write(sample_global_file_string, global_file_write_test_fname):
#     g_init = Global(sample_global_file_string)
#     g_init.write(global_file_write_test_fname)
#     print('global_file_write_test_fname: {}'.format(global_file_write_test_fname))

#     with open(global_file_write_test_fname, 'r') as f:
#         g_final = Global(f)

#     for gp in g_init:
#         assert gp in g_final
