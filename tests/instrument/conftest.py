import pytest
import numina.instrument.assembly as asb

import emirdrp.loader
from emirdrp.testing.create_headers import create_dtu_header_example


@pytest.fixture
def dtu_header_example():
    return create_dtu_header_example()


@pytest.fixture(scope="session")
def component_store():
    """The components of the instrument configurations of EMIR"""
    drp = emirdrp.loader.load_drp()
    return asb.load_paths_store([drp.profiles], [])
