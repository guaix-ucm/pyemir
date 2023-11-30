
import numpy
import math
import pytest
import astropy.units as u
from hypothesis import given

from ..dtu import DtuConf, average, DetectorTranslationUnit
from ..dtuaxis import managed_ndig
from emirdrp.instrument.tests.dtuheader import dtu_fits_header, DTU_HEADER_EXAMPLE


@given(dtu_fits_header())
def test_dtuc(hdr):
    dtuconf = DtuConf.from_header(hdr)
    assert isinstance(dtuconf, DtuConf)


def test_dtu_values():
    dtuconf = DtuConf.from_values(xdtu=1.0, ydtu=1.0, zdtu=1.0)
    assert isinstance(dtuconf, DtuConf)


def test_dtu_values_raises():
    with pytest.raises(ValueError):
        DtuConf.from_values(xdtu=1.0, zdtu=1.0)


@given(dtu_fits_header())
def test_dtu_eq(hdr):
    dtuconf1 = DtuConf.from_header(hdr)
    kwvals = {(k.lower()): v for k, v in hdr.items()}
    dtuconf2 = DtuConf.from_values(**kwvals)
    assert dtuconf1 == dtuconf2


def test_dtuc_shift():
    import emirdrp.instrument.constants as cons

    # Value in microns in DTU coords
    header = DTU_HEADER_EXAMPLE
    expected = [-17.206285211909773, -0.7186057740102463, 0.055023193358977096]
    dtuconf = DtuConf.from_header(header)

    coor_r = dtuconf.coor_r
    assert numpy.allclose(coor_r, expected)

    expected_pix = [0.03992254, -0.95590473,  0.00305684]
    # Value in pixels in image coords
    trans3 = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]  # T3 = T2 * T1

    vec = numpy.dot(trans3, dtuconf.coor_r) * u.micron / cons.EMIR_PIXSIZE
    assert numpy.allclose(vec, expected_pix)


def test_dtur():
    header = DTU_HEADER_EXAMPLE
    dtuconf = DtuConf.from_header(header)
    x_dtu = [-205.679000854492, -24.4878005981445, -463.765991210938]
    x_dtu_r = [-17.206285211909773, -0.7186057740102463, 0.055023193358977096]
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
    header = DTU_HEADER_EXAMPLE
    dtuconf = DtuConf.from_header(header)
    val = dtuconf.describe()
    assert val == expected


@given(dtu_fits_header())
def test_dtu_managed(hdr):
    dtuconf = DtuConf.from_header(hdr)
    with managed_ndig(dtuconf, 4):
        assert dtuconf.get_ndig() == 4

    assert dtuconf.get_ndig() == 3


@given(dtu_fits_header(), dtu_fits_header())
def test_dtu_average(hdr1, hdr2):
    dtuconf1 = DtuConf.from_header(hdr1)
    dtuconf2 = DtuConf.from_header(hdr2)
    dtuconf = average(dtuconf1, dtuconf2)
    assert isinstance(dtuconf, DtuConf)


def test_dtu_average0():
    dtuconf = average()
    assert isinstance(dtuconf, DtuConf)


@given(dtu_fits_header(), dtu_fits_header())
def test_dtu_average2_1(hdr1, hdr2):
    dtuconf1 = DetectorTranslationUnit.from_header(hdr1)
    dtuconf2 = DetectorTranslationUnit.from_header(hdr2)
    dtuconf = average(dtuconf1, dtuconf2)
    assert isinstance(dtuconf, DetectorTranslationUnit)


@given(dtu_fits_header(), dtu_fits_header())
def test_dtu_average2_2(hdr1, hdr2):
    dtuconf1 = DetectorTranslationUnit.from_header(hdr1)
    dtuconf2 = DetectorTranslationUnit.from_header(hdr2)
    dtuconf = average(dtuconf1, dtuconf2)
    for name in ['xaxis', 'yaxis', 'zaxis']:
        calc = 0.5 * (getattr(dtuconf1, name).coor + getattr(dtuconf2, name).coor)
        assert math.isclose(getattr(dtuconf, name).coor, calc)
    for name in ['xaxis', 'yaxis', 'zaxis']:
        calc = 0.5 * (getattr(dtuconf1, name).coor_f + getattr(dtuconf2, name).coor_f)
        assert math.isclose(getattr(dtuconf, name).coor_f, calc)
    for name in ['xaxis', 'yaxis', 'zaxis']:
        calc = 0.5 * (getattr(dtuconf1, name).coor_0 + getattr(dtuconf2, name).coor_0)
        assert math.isclose(getattr(dtuconf, name).coor_0, calc)
