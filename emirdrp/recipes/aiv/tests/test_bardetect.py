
import pytest
import numpy

from ..bardetect import calc_fwhm
from ..bardetect import find_position


@pytest.mark.parametrize("axis, result", [(0,7), (1, 16)])
def test_calc_fwhm(axis, result):

    img = numpy.zeros((100, 100))
    img[42:58,45:52] = 100
    region = (slice(40,60), slice(40, 60))
    fwhm = calc_fwhm(img, region, axis=axis)
    assert fwhm == result


def test_find_position():

    edges = numpy.zeros((100, 100), dtype='int')
    edges[5:12,40] = 1
    edges[5:12,45] = 1

    # Find slit
    pos = find_position(edges, 8.0, 1, 99)
    assert pos == (8, 40.0, 45.0)

    # Dont find it
    pos = find_position(edges, 40, 1, 99)
    assert pos is None