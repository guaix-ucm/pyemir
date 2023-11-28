
import numpy
import pytest

from ..dtuaxis import DtuAxisAdaptor, DTUAxis
from ..dtuconf import apply_on_axis

from emirdrp.instrument.tests.dtuheader import HEADERS


def test_dtuaxis_raise():
    with pytest.raises(ValueError):
        DtuAxisAdaptor("R", 200.0)


def test_dtuaxis_header():

    header = HEADERS[0]

    dtuaxis_x = DtuAxisAdaptor.from_header(header, name='X')

    assert isinstance(dtuaxis_x, DtuAxisAdaptor)

def test_dtuaxis_header2():

    header = HEADERS[0]

    dtuaxis_x = DTUAxis('X')
    dtuaxis_x.configure_with_header(header)

    assert dtuaxis_x.is_configured



@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuaxis_header_raise(hdr):
    with pytest.raises(ValueError):
        DtuAxisAdaptor.from_header(hdr, name='R')


def test_dtu_apply():
    axis1 = DtuAxisAdaptor("X", 100.0)
    axis2 = DtuAxisAdaptor("X", 120.0)

    axis3 = apply_on_axis(numpy.mean, [axis1, axis2])
    assert axis3.coor == 110.0
    assert axis3.coor_f == 1.0
    assert axis3.coor_0 == 0.0
    assert axis3.name == "X"


def test_dtu_apply0():
    """Test on an empty list"""
    axis3 = apply_on_axis(numpy.mean, [])
    assert numpy.isnan(axis3.coor)
    assert numpy.isnan(axis3.coor_f)
    assert numpy.isnan(axis3.coor_0)
    assert axis3.name == "X"
