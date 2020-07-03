
import pytest
import numpy
import astropy.units as u
import astropy.wcs

import emirdrp.instrument.constants as cons
from ..distortions import exvp, pvex
from ..distortions import wcs_exvp, wcs_pvex
from ..distortions import adapt_wcs


def create_wcs(fppa=270 * u.deg):

    import astropy.wcs

    w = astropy.wcs.WCS(naxis=2)
    ipa = cons.EMIR_REF_IPA
    angle = fppa - ipa
    # scale = numpy.cos(numpy.deg2rad(ipa - fppa))
    w.wcs.ctype = ['RA---ZPN', 'DEC--ZPN']
    w.wcs.crval = [86.3712733276122, 28.9182815703162]
    w.wcs.crpix = [1022.68, 1016.70]
    # scale = cons.EMIR_PIXSCALE.to(u.deg / u.pixel).value
    scale1 = 0.194279723 / 3600
    scale2 = 0.194264445 / 3600
    #w.wcs.cdelt = [scale, scale]
    pc11 = scale1 * numpy.cos(angle)
    pc22 = scale2 * numpy.cos(angle)
    pc12 = -scale1 *numpy.sin(angle)
    pc21 = scale2 * numpy.sin(angle)
    w.wcs.cd = [[pc11, pc12],
                [pc21, pc22]]
    w.wcs.set_pv([
        (2, 1, 0.9998), (2, 2, 0.0), (2, 3, 13824.76),
        (2, 4, 0.0), (2, 5, 3491446467)
    ])
    w.wcs.radesys = 'FK5'
    return w


def create_wcs_new():
    """From EMIR image as 2020-07-05"""
    w = astropy.wcs.WCS(naxis=2)
    w.wcs.ctype = ['RA---ZPN', 'DEC--ZPN']
    w.wcs.crval = [297.591888896875, 30.7466139938489]
    w.wcs.crpix = [1022.46378968092, 1016.55073697444]
    pc11 = -5.39644538444764E-05
    pc22 = -5.39644538444764E-05
    pc12 = -5.19905444897492E-08
    pc21 =  5.19905444897492E-08
    w.wcs.cd = [[pc11, pc12],
                [pc21, pc22]]
    w.wcs.set_pv([
        (2, 1, 1), (2, 2, 0.0), (2, 3, 14584.8),
        (2, 4, 6.93534705104727E-310), (2, 5, 1755295542.4)
    ])
    w.wcs.radesys = 'FK5'
    return w


def test_adapt_wcs():
    wcs = create_wcs()
    wcsa = adapt_wcs(wcs, 90.0552, 90.0552)
    assert id(wcsa.wcs) != id(wcs.wcs)


def create_wcs_alt(fppa=270 * u.deg):
    """Candidate alternate WCS to pass to virtual pixels, not working"""
    import astropy.wcs

    w = astropy.wcs.WCS(naxis=2)
    ipa = cons.EMIR_REF_IPA
    angle = fppa - ipa
    # scale = numpy.cos(numpy.deg2rad(ipa - fppa))
    w.wcs.ctype = ['', '']
    w.wcs.crval = [86.3712733276122, 28.9182815703162]
    w.wcs.crpix = [1022.68, 1016.70]
    # scale = cons.EMIR_PIXSCALE.to(u.deg / u.pixel).value
    scale1 = 0.194279723 / 3600
    scale2 = 0.194264445 / 3600
    #w.wcs.cdelt = [scale, scale]
    pc11 = scale1 * numpy.cos(angle)
    pc22 = scale2 * numpy.cos(angle)
    pc12 = -scale1 *numpy.sin(angle)
    pc21 = scale2 * numpy.sin(angle)
    w.wcs.cd = [[pc11, pc12],
                [pc21, pc22]]
    w.wcs.radesys = 'FK5'
    return w


def test_distortions_ex1():

    e_x1 = [-14.95369787,    88.60095332,   999.76814675,  2019.84941119]
    e_y1 = [1202.73558767,   898.46492016,   192.1974284,    -26.51887593]
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
    p_x1 = [16.20132521,   110.97939823,  1000.22549019,  1962.80728145]
    p_y1 = [1197.39342201,   901.47856688,   207.58843519,    33.71359561]
    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]

    x1, y1 = pvex(x0, y0)

    assert numpy.allclose(x1, p_x1)
    assert numpy.allclose(y1, p_y1)


def test_distortions_pc2():
    """Inverse transformation, not very precise"""

    e_x1 = [-14.95369787,    88.60095332,   999.76814675,  2019.84941119]
    e_y1 = [1202.73558767,   898.46492016,   192.1974284,    -26.51887593]

    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]
    # x0 = [0.992795185,   100.002590,   1000.00025,   1989.91083]
    # y0 = [1200.00123542, 900.00034883, 200.0082549, 5.09416157]

    p_x0, p_y0 = pvex(e_x1, e_y1)

    assert numpy.allclose(x0, p_x0, rtol=2e-2)
    assert numpy.allclose(y0, p_y0, rtol=2e-2)


def test_distortions_wcs_ex():

    w = create_wcs()

    e_x1 = [ -15.73469527,   88.61357533, 1000.58521175, 2024.25873613]
    e_y1 = [1201.96995355,  897.6454293 ,  192.4567596 ,  -28.81253997]
    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]
    x1, y1 = wcs_exvp(w, x0, y0)

    assert numpy.allclose(x1, e_x1)
    assert numpy.allclose(y1, e_y1)


def test_distortions_wcs_pc():

    w = create_wcs()

    p_x1 = [  16.87732924,  110.93407453,  999.42297684, 1959.34745267]
    p_y1 = [1198.15178941,  902.27528867,  207.32346469,   35.16785173]
    x0 = [1, 100, 1000, 1990]
    y0 = [1200, 900, 200, 5]

    x1, y1 = wcs_pvex(w, x0, y0)

    assert numpy.allclose(x1, p_x1)
    assert numpy.allclose(y1, p_y1)


def test_consistency_wcs():
    x0 = [1, 1, 2000, 2000]
    y0 = [1, 2000, 2000, 1]

    w = create_wcs()

    v_x1 = [  33.61421368,   34.25887803, 1971.56801359, 1968.54383359]
    v_y1 = [  35.31658127, 1969.82615199, 1969.55586265,   31.79447425]

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

    v_x1 = [  33.61421368,   34.25887803, 1971.56801359, 1968.54383359]
    v_y1 = [  35.31658127, 1969.82615199, 1969.55586265,   31.79447425]

    ra_i, dec_i = wcs1.all_pix2world(x0, y0, 1)
    x1, y1 = wcs2.all_world2pix(ra_i, dec_i, 1)

    assert numpy.allclose(x1, v_x1)
    assert numpy.allclose(y1, v_y1)
