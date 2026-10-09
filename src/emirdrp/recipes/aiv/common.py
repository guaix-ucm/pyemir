#
# Copyright 2013-2024 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


"""AIV Recipes for EMIR"""

import logging

import numpy
from scipy import ndimage


from numina.array.fwhm import compute_fwhm_2d_simple
from numina.array.utils import expand_region

_logger = logging.getLogger(__name__)


# returns y,x
# Background in an annulus, mode is HSM
# returns x,y
def shape_of_slices(tup_of_s):
    return tuple(m.stop - m.start for m in tup_of_s)


def normalize(data):
    b = data.max()
    a = data.min()

    if b != a:
        data_22 = (2 * data - b - a) / (b - a)
    else:
        data_22 = data - b
    return data_22


def normalize_raw(arr):
    """Rescale float image between 0 and 1

    This is an extension of image_as_float of scikit-image
    when the original image was uint16 but later was
    processed and transformed to float32

    Parameters
    ----------
    arr; ndarray

    Returns
    -------
      A ndarray mapped between 0 and 1
    """

    # FIXME: use other limits acording to original arr.dtype
    # This applies only to uint16 images
    # As images were positive, the range is 0,1

    return numpy.clip(arr / 65535.0, 0.0, 1.0)


def char_slit(data, regions, box_increase=3, slit_size_ratio=4.0):

    result = []

    for r in regions:
        _logger.debug("initial region %s", r)
        oshape = shape_of_slices(r)

        ratio = oshape[0] / oshape[1]
        if (slit_size_ratio > 0) and (ratio < slit_size_ratio):
            _logger.debug("this is not a slit, ratio=%f", ratio)
            continue

        _logger.debug("initial shape %s", oshape)
        _logger.debug("ratio %f", ratio)
        rp = expand_region(r, box_increase, box_increase, start=0, stop=2048)
        _logger.debug("expanded region %r", rp)
        ref = rp[0].start, rp[1].start
        _logger.debug("reference point %r", ref)

        datas = data[rp]

        c = ndimage.center_of_mass(datas)

        fc = datas.shape[0] // 2
        cc = datas.shape[1] // 2
        _logger.debug("%d %d %d %d", fc, cc, c[0], c[1])

        _peak, fwhm_x, fwhm_y = compute_fwhm_2d_simple(datas, c[1], c[0])

        _logger.debug("x=%f y=%f", c[1] + ref[1], c[0] + ref[0])
        _logger.debug("fwhm_x %f fwhm_y %f", fwhm_x, fwhm_y)

        # colrow = ref[1] + cc + 1, ref[0] + fc + 1

        result.append([c[1] + ref[1] + 1, c[0] + ref[0] + 1, fwhm_x, fwhm_y])

    return result
