
import pytest

import numpy.random
import scipy.ndimage
import numina.array.utils as utils

from ..corr import offsets_from_crosscor


@pytest.fixture(scope='module')
def images():
    numpy.random.seed(1202929)
    result = []
    shape = (1000, 1000)
    off = [(500, 500), (480, 490), (505, 500), (500, 500)]

    for idx, o in enumerate(off):
        data = numpy.zeros(shape) + 1000
        data[o] = 50000
        data = scipy.ndimage.gaussian_filter(data, 2.0)
        data = numpy.random.normal(data, 30.0)
        result.append(data)

    return result


@pytest.mark.parametrize("refine", [True, False])
@pytest.mark.parametrize("order", ['ij', 'xy'])
def test_coor(images, order, refine):

    arrs = images

    shape = arrs[0].shape
    xref_cross = shape[1] // 2
    yref_cross = shape[0] // 2
    box = 50

    if refine:
        result = numpy.array([
            [0.00000000e+00, 0.00000000e+00],
            [-1.99935691e+01, -1.00047803e+01],
            [4.98913940e+00, -1.41737914e-02],
            [2.93443492e-02, -3.53881751e-02]
        ]
        )
    else:
        result = numpy.array([
            [0.0, 0.0],
            [-20.0, -10.0],
            [5.0, 0.0],
            [0.0, 0.0]
        ]
        )

    if order == 'xy':
        result = result[:, ::-1]

    region = utils.image_box2d(xref_cross, yref_cross, shape, (box, box))
    computed = offsets_from_crosscor(arrs, region, refine=refine, order=order)

    assert numpy.allclose(result, computed)


def test_coor_raises(images):

    arrs = images

    shape = arrs[0].shape
    xref_cross = shape[1] // 2
    yref_cross = shape[0] // 2
    box = 50

    region = utils.image_box2d(xref_cross, yref_cross, shape, (box, box))
    with pytest.raises(ValueError):
        offsets_from_crosscor(arrs, region, order="sksjd")
