import numpy
from astropy.modeling.functional_models import Gaussian2D


def create_scene_noise(size, background=0, std=100, pos=3, peak=0, dtype=None):
    data = numpy.random.normal(background, std, size)
    if pos is not None:
        data[pos] += peak
    return data.astype(dtype)


def create_scene_23(size, pos=3, peak=0, dtype=None):
    data = numpy.zeros(size, dtype=dtype)
    if pos is not None:
        data[pos] += peak
    return data.astype(dtype)


def create_scene_4847(val=0):
    data = numpy.zeros((10, 10)) + val
    return data


def create_scene_4829(size, wcs_obj, background=0, std=100.0, dtype=None):
    data = numpy.random.normal(background, std, size)
    y, x = numpy.mgrid[: size[1], : size[0]]
    model = Gaussian2D(
        amplitude=30 * background,
        x_mean=wcs_obj.wcs.crpix[0] - 3.5,
        y_mean=wcs_obj.wcs.crpix[1] - 12.8,
        x_stddev=3.0,
        y_stddev=4.0,
    )
    data += model(x, y)
    return data.astype(dtype)
