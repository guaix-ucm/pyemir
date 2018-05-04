
import numpy
from emirdrp.instrument.distortions import exvp, pvex


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
