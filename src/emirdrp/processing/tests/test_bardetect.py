
import pytest
import numpy

from emirdrp.processing.bardetect import calc_fwhm
from emirdrp.processing.bardetect import find_position

@pytest.mark.parametrize("axis, result", [(0,7), (1, 16)])
def test_calc_fwhm(axis, result):

    img = numpy.zeros((100, 100))
    img[42:58,45:52] = 100
    region = (slice(40,60), slice(40, 60))
    fwhm = calc_fwhm(img, region, axis=axis)
    assert fwhm == result


def test_find_position_2():
    """Two bars"""

    edges = numpy.zeros((20, 20), dtype='bool')

    edges[3:5,5] = True
    edges[5:8,5] = True

    edges[3:5,12] = True
    edges[5:8,12] = True

    m = find_position(edges, 5, 2, 19, total=5)
    m2 = sorted(m, key=lambda cen: cen[0])

    assert len(m2) == 2

    assert m2[0] == (5.0, 5.0)
    assert m2[1] == (12.0, 5.0)


def test_find_position_0():
    """No edges"""
    edges = numpy.zeros((20, 20), dtype='bool')

    m = find_position(edges, 5, 2, 19, total=5)

    assert len(m) == 0


def test_find_position_1():
    """Only 1 edge."""
    edges = numpy.zeros((20, 20), dtype='bool')

    edges[3:5,5] = True
    edges[5:8,4] = True

    m = find_position(edges, 5, 2, 19, total=5)

    res = (4.4000000000000004, 5.0)

    assert len(m) == 1

    assert m[0] == res
