import pytest

from emirdrp.testing.create_headers import create_dtu_header_example


@pytest.fixture
def dtu_header_example():
    return create_dtu_header_example()
