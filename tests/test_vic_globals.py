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

def test_global_init():
    g = Global('/tmp/glb_peyto_base_save_state.txt')
    print(str(g))
