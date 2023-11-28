
import numpy
import pytest
import hypothesis.strategies as st
from hypothesis import given

from ..dtuaxis import DtuAxisAdaptor, DTUAxis
from ..dtuconf import apply_on_axis

from emirdrp.instrument.tests.dtuheader import dtu_fits_header


def test_dtuaxis_raise():
    with pytest.raises(ValueError):
        DtuAxisAdaptor("R", 200.0)


@given(dtu_fits_header())
def test_dtuaxis_header(hdr):
    dtuaxis_x = DtuAxisAdaptor.from_header(hdr, name='X')
    assert isinstance(dtuaxis_x, DtuAxisAdaptor)


@given(dtu_fits_header(), st.sampled_from(['X', 'Y', 'Z']))
def test_dtuaxis_header2(hdr, name):
    dtuaxis = DTUAxis(name)
    dtuaxis.configure_with_header(hdr)
    assert dtuaxis.is_configured


@given(dtu_fits_header())
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
