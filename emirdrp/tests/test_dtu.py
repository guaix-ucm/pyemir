
import numpy
import pytest

from emirdrp.datamodel import DtuConf, DtuAxis


# header fragment

def dtu_header():
    header = {}
    header['XDTU_F'] = 0.922918
    header['YDTU_F'] = 0.971815
    header['ZDTU_F'] = 1.0
    header['XDTU'] = -205.679000854492
    header['YDTU'] = -24.4878005981445
    header['ZDTU'] = -463.765991210938
    header['XDTU_0'] = -205.651000976562
    header['YDTU_0'] = -24.4794006347656
    header['ZDTU_0'] = -463.821014404297
    return header


def test_dtuc():

    header = dtu_header()

    dtuconf = DtuConf.from_header(header)

    assert isinstance(dtuconf, DtuConf)


def test_dtuc_shift():
    from emirdrp.core import EMIR_PIXSCALE

    # Value in microns in DTU coords
    header = dtu_header()
    expected = [-0.030338424356226815, -0.008643582758960445, 0.055023193358977096]
    dtuconf = DtuConf.from_header(header)

    coor_r = dtuconf.coor_r
    assert numpy.allclose(coor_r, expected)

    expected_pix = [0.0004802, -0.00168547, 0.00305684]
    # Value in pixels in image coords
    trans3 = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]  # T3 = T2 * T1

    vec = numpy.dot(trans3, dtuconf.coor_r) / EMIR_PIXSCALE
    assert numpy.allclose(vec, expected_pix)


def test_dtuaxis_raise():
    with pytest.raises(ValueError):
        DtuAxis("R", 200.0)


def test_dtuaxis_header():

    header = dtu_header()

    dtuaxis_x = DtuAxis.from_header(header, name='X')

    assert isinstance(dtuaxis_x, DtuAxis)


def test_dtuaxis_header_raise():

    header = dtu_header()

    with pytest.raises(ValueError):
        DtuAxis.from_header(header, name='R')
