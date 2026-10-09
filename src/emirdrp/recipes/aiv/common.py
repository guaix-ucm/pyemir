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
import math

import numpy
import scipy.interpolate as itpl
import scipy.optimize as opz
from scipy import ndimage

from astropy.modeling import models, fitting

from numina.array.recenter import centering_centroid
from numina.array.utils import image_box
from numina.array.fwhm import compute_fwhm_2d_simple
from numina.array.utils import expand_region

_logger = logging.getLogger(__name__)


# returns y,x
def compute_fwhm(img, center):
    X = numpy.arange(img.shape[0])
    Y = numpy.arange(img.shape[1])

    bb = itpl.RectBivariateSpline(X, Y, img)
    # We assume that the peak is in the center...
    peak = bb(*center)[0, 0]

    def f1(x):
        return bb(x, center[1]) - 0.5 * peak

    def f2(y):
        return bb(center[0], y) - 0.5 * peak

    def compute_fwhm_1(U, V, fun, center):

        cp = int(math.floor(center + 0.5))

        # Min on the rigth
        r_idx = V[cp:].argmin()
        u_r = U[cp + r_idx]

        if V[cp + r_idx] > 0.5 * peak:
            # FIXME: we have a problem
            # brentq will raise anyway
            pass

        sol_r = opz.brentq(fun, center, u_r)

        # Min in the left
        rV = V[cp - 1 :: -1]
        rU = U[cp - 1 :: -1]
        l_idx = rV.argmin()
        u_l = rU[l_idx]
        if rV[l_idx] > 0.5 * peak:
            # FIXME: we have a problem
            # brentq will raise anyway
            pass

        sol_l = opz.brentq(fun, u_l, center)
        fwhm = sol_r - sol_l
        return fwhm

    U = X
    V = bb.ev(U, [center[1] for _ in U])
    fwhm_x = compute_fwhm_1(U, V, f1, center[0])

    U = Y
    V = bb.ev([center[0] for _ in U], U)

    fwhm_y = compute_fwhm_1(U, V, f2, center[1])

    return center[0], center[1], peak, fwhm_x, fwhm_y


# Background in an annulus, mode is HSM
def compute_fwhm_global(data, center, box):
    sl = image_box(center, data.shape, box)
    raster = data[sl]

    background = raster.min()
    braster = raster - background

    newc = center[0] - sl[0].start, center[1] - sl[1].start
    try:
        res = compute_fwhm(braster, newc)
        return (res[1] + sl[1].start, res[0] + sl[0].start, res[2], res[3], res[4])
    except ValueError as error:
        _logger.warning("%s", error)
        return center[1], center[0], -99.0, -99.0, -99.0
    except Exception as error:
        _logger.warning("%s", error)
        return center[1], center[0], -199.0, -199.0, -199.0


# returns x,y
def gauss_model(data, center_r):
    sl = image_box(center_r, data.shape, box=(4, 4))
    raster = data[sl]

    # background
    background = raster.min()

    b_raster = raster - background

    new_c = center_r[0] - sl[0].start, center_r[1] - sl[1].start

    yi, xi = numpy.indices(b_raster.shape)

    g = models.Gaussian2D(
        amplitude=b_raster.max(),
        x_mean=new_c[1],
        y_mean=new_c[0],
        x_stddev=1.0,
        y_stddev=1.0,
    )
    f1 = fitting.LevMarLSQFitter()  # @UndefinedVariable
    t = f1(g, xi, yi, b_raster)

    mm = (
        t.x_mean.value + sl[1].start,
        t.y_mean.value + sl[0].start,
        t.amplitude.value,
        t.x_stddev.value,
        t.y_stddev.value,
    )
    return mm


def recenter_char(data, centers_i, recenter_maxdist, recenter_nloop, recenter_half_box, do_recenter):

    # recentered values
    centers_r = numpy.empty_like(centers_i)
    # Ignore certain pinholes
    compute_mask = numpy.ones((centers_i.shape[0],), dtype="bool")
    status_array = numpy.ones((centers_i.shape[0],), dtype="int")

    for idx, (xi, yi) in enumerate(centers_i):
        # A failsafe
        _logger.info("for pinhole %i", idx)
        _logger.info("center is x=%7.2f y=%7.2f", xi, yi)
        if xi > data.shape[1] - 5 or xi < 5 or yi > data.shape[0] - 5 or yi < 5:
            _logger.info("pinhole too near to the border")
            compute_mask[idx] = False
            centers_r[idx] = xi, yi
            status_array[idx] = 0
        else:
            if do_recenter and (recenter_maxdist > 0.0):
                _ = centering_centroid(
                    data,
                    xi,
                    yi,
                    box=recenter_half_box,
                    maxdist=recenter_maxdist,
                    nloop=recenter_nloop,
                )
                xc, yc, _back, status, msg = centering_centroid(
                    data,
                    xi,
                    yi,
                    box=recenter_half_box,
                    maxdist=recenter_maxdist,
                    nloop=recenter_nloop,
                )
                _logger.info("new center is x=%7.2f y=%7.2f", xc, yc)
                # Log in X,Y format
                _logger.debug("recenter message: %s", msg)
                centers_r[idx] = xc, yc
                status_array[idx] = status
            else:
                centers_r[idx] = xi, yi
                status_array[idx] = 0

    return centers_r, compute_mask, status_array


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
