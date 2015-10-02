
import pytest
import numpy

from ..bardetect import calc_fwhm

@pytest.mark.parametrize("axis, result", [(0,7), (1, 16)])
def test_calc_fwhm(axis, result):

    img = numpy.zeros((100, 100))
    img[42:58,45:52] = 100
    region = (slice(40,60), slice(40, 60))
    fwhm = calc_fwhm(img, region, axis=axis)
    assert fwhm == result

