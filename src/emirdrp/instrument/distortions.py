#
# Copyright 2008-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numpy
import astropy.units as U

import emirdrp.instrument.constants as cons

EMIR_PLATESCALE_RADS = cons.EMIR_PIXSCALE.to(U.rad / U.pixel).value


def adapt_wcs(wcs, rotang_deg, ipa_deg):
    """Adapt WCS to use in distortions"""
    wcsa = wcs.deepcopy()
    wcsa.wcs.crval = [0, 0]
    cdelt = numpy.sqrt(numpy.linalg.det(wcs.wcs.cd))
    angle = numpy.arctan2(wcs.wcs.cd[1, 0], wcs.wcs.cd[0, 0])
    nangle = angle + (ipa_deg - rotang_deg) / 180.0 * numpy.pi
    pc11 = cdelt * numpy.cos(nangle)
    pc22 = cdelt * numpy.cos(nangle)
    pc12 = -cdelt * numpy.sin(nangle)
    pc21 = cdelt * numpy.sin(nangle)
    wcsa.wcs.cd = [[pc11, pc12],
                   [pc21, pc22]]
    return wcsa


def exvp(pos_x, pos_y):
    """Convert virtual pixel to real pixel

    This function is deprecated, use wcs_exvp instead
    """
    pos_x = numpy.asarray(pos_x)
    pos_y = numpy.asarray(pos_y)
    # convert virtual pixel to real pixel
    # convert world coordinate to pixel
    center = [1024.5, 1024.5]
    cf = EMIR_PLATESCALE_RADS
    pos_base_x = pos_x - center[0]
    pos_base_y = pos_y - center[1]
    ra = numpy.hypot(pos_base_x, pos_base_y)
    thet = numpy.arctan2(pos_base_y, pos_base_x)

    r  = cf * ra

    rr1 = 1 + 14606.7 * r**2 + 1739716115.1 * r**4

    nx1 = rr1 * ra * numpy.cos(thet) + center[0]
    ny1 = rr1 * ra * numpy.sin(thet) + center[1]
    return nx1, ny1


def wcs_exvp(wcs, xv, yv):
    """Convert virtual pixel to real pixel"""

    # These factors are equivalent to 1 / CDELT1 1 / CDELT2
    cdelt1 = numpy.hypot(wcs.wcs.cd[0, 0], wcs.wcs.cd[0, 1])
    cdelt2 = numpy.hypot(wcs.wcs.cd[1, 0], wcs.wcs.cd[1, 1])
    #
    ra_i = (wcs.wcs.crpix[0] - xv) * cdelt1
    dec_i = (wcs.wcs.crpix[1] - yv) * cdelt2
    # keep old value of crval
    oldcrval = wcs.wcs.crval.copy()
    wcs.wcs.crval = [0, 0]
    xr, yr = wcs.all_world2pix(ra_i, dec_i, 1)
    # reset crval
    wcs.wcs.crval = oldcrval
    return xr, yr


def pvex(pos_x, pos_y):
    """Convert real pixel to virtual pixel


    This function is deprecated, use wcs_pvex instead
    """
    pos_x = numpy.asarray(pos_x)
    pos_y = numpy.asarray(pos_y)

    center = [1024.5, 1024.5]
    cf = EMIR_PLATESCALE_RADS
    pos_base_x = pos_x - center[0]
    pos_base_y = pos_y - center[1]
    ra = numpy.hypot(pos_base_x, pos_base_y)
    thet = numpy.arctan2(pos_base_y, pos_base_x)

    r  = cf * ra
    rr1 = 1.000051 - 14892 * r**2 - 696254464 * r**4

    nx1 = rr1 * ra * numpy.cos(thet) + center[0]
    ny1 = rr1 * ra * numpy.sin(thet) + center[1]
    return nx1, ny1


def wcs_pvex(wcs, xr, yr):
    """Convert real pixel to virtual pixel"""

    # These factors are equivalent to 1 / CDELT1 1 / CDELT2
    cdelt1 = numpy.hypot(wcs.wcs.cd[0, 0], wcs.wcs.cd[0, 1])
    cdelt2 = numpy.hypot(wcs.wcs.cd[1, 0], wcs.wcs.cd[1, 1])
    cfax = 1 / cdelt1
    cfay = 1 / cdelt2

    # keep old value of crval
    oldcrval = wcs.wcs.crval.copy()
    wcs.wcs.crval = [0, 0]
    ra_i, dec_i = wcs.all_pix2world(xr, yr, 1)
    # wrap angles around 180
    # [-180, 180)
    ra_i = numpy.mod(ra_i - 180, 360.0) - 180

    xv = wcs.wcs.crpix[0] - ra_i * cfax
    yv = wcs.wcs.crpix[1] - dec_i * cfay
    # restore old value of crval
    wcs.wcs.crval = oldcrval
    return xv, yv


def pix2virt(pos, origin=1):
    """

    This function is deprecated, use wcs_pix2virt instead
    """
    ddef_o = 1
    off = ddef_o - origin
    pos = numpy.atleast_2d(pos) + off
    nx, ny = pvex(pos[:, 0], pos[:, 1])
    res = numpy.stack((nx, ny), axis=1)
    return res - off


def wcs_pix2virt(wcs, pos, origin=1):
    ddef_o = 1
    off = ddef_o - origin
    pos = numpy.atleast_2d(pos) + off
    nx, ny = wcs_pvex(wcs, pos[:, 0], pos[:, 1])
    res = numpy.stack((nx, ny), axis=1)
    return res - off
