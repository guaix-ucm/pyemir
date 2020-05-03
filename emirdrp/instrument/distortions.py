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


def exvp(pos_x, pos_y):
    """Convert virtual pixel to real pixel"""
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


def pvex(pos_x, pos_y):
    """Convert real pixel to virtual pixel"""
    # convert real pixel to virtual pixel.
    # convert pixel to world coordinate
    pos_x = numpy.asarray(pos_x)
    pos_y = numpy.asarray(pos_y)

    center = [1024.5, 1024.5]
    cf  = EMIR_PLATESCALE_RADS
    pos_base_x = pos_x - center[0]
    pos_base_y = pos_y - center[1]
    ra = numpy.hypot(pos_base_x, pos_base_y)
    thet = numpy.arctan2(pos_base_y, pos_base_x)

    r  = cf * ra
    rr1 = 1.000051 - 14892 * r**2 - 696254464 * r**4

    nx1 = rr1 * ra * numpy.cos(thet) + center[0]
    ny1 = rr1 * ra * numpy.sin(thet) + center[1]
    return nx1, ny1


def pix2virt(pos, origin=1):
    ddef_o = 1
    off = ddef_o - origin
    pos = numpy.atleast_2d(pos) + off
    nx, ny = pvex(pos[:, 0], pos[:, 1])
    res = numpy.stack((nx, ny), axis=1)
    return res - off
