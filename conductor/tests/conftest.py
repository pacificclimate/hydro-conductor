import io
from pkg_resources import resource_stream

import pytest

@pytest.fixture(scope="module")
def sample_global_file_string():
    stream = resource_stream('conductor', 'tests/input/global.txt')
    return io.TextIOWrapper(stream)
