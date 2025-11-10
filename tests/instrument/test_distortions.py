import pytest
import numpy

from emirdrp.instrument.distortions import exvp, pvex
from emirdrp.instrument.distortions import wcs_exvp, wcs_pvex
from emirdrp.instrument.distortions import adapt_wcs
from emirdrp.testing.create_wcs import create_wcs, create_wcs_alt


def test_adapt_wcs():
    wcs = create_wcs()
    wcsa = adapt_wcs(wcs, 90.0552, 90.0552)
    assert id(wcsa.wcs) != id(wcs.wcs)


def test_distortions_ex1():

    e_x1 = [-14.95369787, 88.60095332, 999.76814675, 2019.84941119]
    e_y1 = [1202.73558767, 898.46492016, 192.1974284, -26.51887593]
    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]
    x1, y1 = exvp(x0, y0)

    assert numpy.allclose(x1, e_x1)
    assert numpy.allclose(y1, e_y1)


def test_distortions_ex2():
    """Inverse transformation, not very precise"""

    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]
    p_x1 = [16.20132521, 110.97939823, 1000.22549019, 1962.80728145]
    p_y1 = [1197.39342201, 901.47856688, 207.58843519, 33.71359561]
    x1, y1 = exvp(p_x1, p_y1)

    assert numpy.allclose(x1, x0, rtol=2e-2)
    assert numpy.allclose(y1, y0, rtol=2e-2)


def test_distortions_pc1():
    p_x1 = [16.20132521, 110.97939823, 1000.22549019, 1962.80728145]
    p_y1 = [1197.39342201, 901.47856688, 207.58843519, 33.71359561]
    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]

    x1, y1 = pvex(x0, y0)

    assert numpy.allclose(x1, p_x1)
    assert numpy.allclose(y1, p_y1)


def test_distortions_pc2():
    """Inverse transformation, not very precise"""

    e_x1 = [-14.95369787, 88.60095332, 999.76814675, 2019.84941119]
    e_y1 = [1202.73558767, 898.46492016, 192.1974284, -26.51887593]

    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]
    # x0 = [0.992795185,   100.002590,   1000.00025,   1989.91083]
    # y0 = [1200.00123542, 900.00034883, 200.0082549, 5.09416157]

    p_x0, p_y0 = pvex(e_x1, e_y1)

    assert numpy.allclose(x0, p_x0, rtol=2e-2)
    assert numpy.allclose(y0, p_y0, rtol=2e-2)


def test_distortions_wcs_ex():

    w = create_wcs()

    e_x1 = [-15.73469527, 88.61357533, 1000.58521175, 2024.25873613]
    e_y1 = [1201.96995355, 897.6454293, 192.4567596, -28.81253997]
    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]
    x1, y1 = wcs_exvp(w, x0, y0)

    assert numpy.allclose(x1, e_x1)
    assert numpy.allclose(y1, e_y1)


def test_distortions_wcs_pc():

    w = create_wcs()

    p_x1 = [16.87732924, 110.93407453, 999.42297684, 1959.34745267]
    p_y1 = [1198.15178941, 902.27528867, 207.32346469, 35.16785173]
    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]

    x1, y1 = wcs_pvex(w, x0, y0)

    assert numpy.allclose(x1, p_x1)
    assert numpy.allclose(y1, p_y1)


def test_consistency_wcs():
    x0 = [1, 1, 2000, 2000]
    y0 = [1, 2000, 2000, 1]

    w = create_wcs()

    v_x1 = [33.61421368, 34.25887803, 1971.56801359, 1968.54383359]
    v_y1 = [35.31658127, 1969.82615199, 1969.55586265, 31.79447425]

    x1, y1 = wcs_pvex(w, x0, y0)

    x2, y2 = wcs_exvp(w, x1, y1)

    assert numpy.allclose(x1, v_x1)
    assert numpy.allclose(y1, v_y1)

    assert numpy.allclose(x2, x0)
    assert numpy.allclose(y2, y0)


@pytest.mark.xfail(reason="feature not yet working")
def test_consistency2_wcs():
    """check converting Real<->Virtual using 2 WCS structures"""
    x0 = [1, 1, 2000, 2000]
    y0 = [1, 2000, 2000, 1]

    wcs1 = create_wcs()
    wcs2 = create_wcs_alt()

    v_x1 = [33.61421368, 34.25887803, 1971.56801359, 1968.54383359]
    v_y1 = [35.31658127, 1969.82615199, 1969.55586265, 31.79447425]

    ra_i, dec_i = wcs1.all_pix2world(x0, y0, 1)
    x1, y1 = wcs2.all_world2pix(ra_i, dec_i, 1)

    assert numpy.allclose(x1, v_x1)
    assert numpy.allclose(y1, v_y1)
