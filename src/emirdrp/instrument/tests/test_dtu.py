
import numpy
import pytest
import astropy.units as U

from ..dtuconf import DtuConf, DtuAxis
from ..dtuconf import managed_ndig
from ..dtuconf import average, apply_on_axis
from .dtuheader import HEADERS


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuc(hdr):
    dtuconf = DtuConf.from_header(hdr)
    assert isinstance(dtuconf, DtuConf)


def test_dtu_values():
    dtuconf = DtuConf.from_values(xdtu=1.0, ydtu=1.0, zdtu=1.0)
    assert isinstance(dtuconf, DtuConf)


def test_dtu_values_raises():
    with pytest.raises(ValueError):
        DtuConf.from_values(xdtu=1.0, zdtu=1.0)


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtu_eq(hdr):
    dtuconf1 = DtuConf.from_header(hdr)
    kwvals = {(k.lower()): v for k, v in hdr.items()}
    dtuconf2 = DtuConf.from_values(**kwvals)
    assert dtuconf1 == dtuconf2


def test_dtuc_shift():
    import emirdrp.instrument.constants as cons

    # Value in microns in DTU coords
    header = HEADERS[0]
    expected = [-17.206285211909773, -0.7186057740102463, 0.055023193358977096]
    dtuconf = DtuConf.from_header(header)

    coor_r = dtuconf.coor_r
    assert numpy.allclose(coor_r, expected)

    expected_pix = [ 0.03992254, -0.95590473,  0.00305684]
    # Value in pixels in image coords
    trans3 = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]  # T3 = T2 * T1

    vec = numpy.dot(trans3, dtuconf.coor_r) * U.micron / cons.EMIR_PIXSIZE
    assert numpy.allclose(vec, expected_pix)


def test_dtuaxis_raise():
    with pytest.raises(ValueError):
        DtuAxis("R", 200.0)


def test_dtuaxis_header():

    header = HEADERS[0]

    dtuaxis_x = DtuAxis.from_header(header, name='X')

    assert isinstance(dtuaxis_x, DtuAxis)


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuaxis_header_raise(hdr):
    with pytest.raises(ValueError):
        DtuAxis.from_header(hdr, name='R')


def test_dtur():
    header = HEADERS[0]
    dtuconf = DtuConf.from_header(header)
    x_dtu = [-205.679000854492, -24.4878005981445, -463.765991210938]
    x_dtu_r =  [-17.206285211909773, -0.7186057740102463, 0.055023193358977096]

    assert numpy.allclose(x_dtu, dtuconf.coor)

    assert numpy.allclose(x_dtu_r, dtuconf.coor_r)


def test_dtu_formatter():
    expected = (
        "<DtuConf instance>\n"
        "- XDTU..: -205.679\n"
        "- YDTU..:  -24.488\n"
        "- ZDTU..: -463.766\n"
        "- XDTU_0: -205.651\n"
        "- YDTU_0:  -24.479\n"
        "- ZDTU_0: -463.821\n"
        "- XDTU_F:    0.923\n"
        "- YDTU_F:    0.972\n"
        "- ZDTU_F:    1.000\n"
    )
    header = HEADERS[0]
    dtuconf = DtuConf.from_header(header)
    val = dtuconf.describe()
    assert val == expected


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtu_managed(hdr):
    dtuconf = DtuConf.from_header(hdr)
    with managed_ndig(dtuconf, 4):
        assert dtuconf.get_ndig() == 4

    assert dtuconf.get_ndig() == 3


def test_dtu_apply():
    axis1 = DtuAxis("X", 100.0)
    axis2 = DtuAxis("X", 120.0)

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


def test_dtu_average():
    header1 = HEADERS[0]
    dtuconf1 = DtuConf.from_header(header1)
    header2 = HEADERS[1]
    dtuconf2 = DtuConf.from_header(header2)

    dtuconf = average(dtuconf1, dtuconf2)

    assert isinstance(dtuconf, DtuConf)


def test_dtu_average0():
    dtuconf = average()
    assert isinstance(dtuconf, DtuConf)
