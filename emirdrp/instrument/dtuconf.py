#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""DTU configuration class"""

from __future__ import division

import contextlib

import numpy
import astropy.wcs


@contextlib.contextmanager
def managed_ndig(obj, ndig):
    """Context manager to handle ndig"""
    old_ndig = obj.get_ndig()
    obj.set_ndig(ndig)
    try:
        yield ndig
    finally:
        # Code to release resource, e.g.:
        obj.set_ndig(old_ndig)


class DtuAxis(object):
    """Represents one DTU axis of movement"""
    def __init__(self, name, coor, coor_f=1.0, coor_0=0.0):
        if name.lower() not in ['x', 'y', 'z']:
            raise ValueError('"name" must be "X", "Y" or "Z')
        self.name = name
        self.coor = coor
        self.coor_f = coor_f
        self.coor_0 = coor_0
        self.ndig = 3

    @property
    def coor_r(self):
        return (self.coor / self.coor_f) - self.coor_0

    @classmethod
    def from_header(cls, hdr, name):
        if name.lower() not in ['x', 'y', 'z']:
            raise ValueError('"name" must be "X", "Y" or "Z')
        coor = hdr['{}DTU'.format(name)]
        coor_f = hdr.get('{}DTU_F'.format(name), 1.0)
        coor_0 = hdr.get('{}DTU_0'.format(name), 0.0)
        return DtuAxis(name, coor, coor_f, coor_0)

    def allclose(self, other, rtol=1e-05, atol=1e-08, equal_nan=False):
        import numpy
        a = [self.coor, self.coor_f, self.coor_0]
        b = [other.coor, other.coor_f, other.coor_0]
        return numpy.allclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)

    def closeto(self, other, abserror):
        return self.allclose(other, rtol=0.0, atol=abserror, equal_nan=False)

    def set_ndig(self, ndig):
        self.ndig = ndig

    def get_ndig(self):
        return self.ndig

    def __eq__(self, other):
        # note: set the precision (number of decimal places) to the same
        # number employed in __str__() function above to print out member
        # values
        # FIXME: this is eqv to allclose with atol=1e-4?
        if isinstance(other, DtuAxis):
            result = \
                self.name.lower() == other.name.lower() and \
                (round(self.coor, self.ndig) == round(other.coor, self.ndig)) and \
                (round(self.coor_f, self.ndig) == round(other.coor_f, self.ndig)) and \
                (round(self.coor_0, self.ndig) == round(other.coor_0, self.ndig))

            return result
        return NotImplemented

    def __ne__(self, other):
        return not self == other


class DtuConf(object):
    """Represents the three DTU axes of movement"""
    def __init__(self, xaxis, yaxis, zaxis):
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.zaxis = zaxis

    def set_ndig(self, ndig):
        self.xaxis.set_ndig(ndig)
        self.yaxis.set_ndig(ndig)
        self.zaxis.set_ndig(ndig)

    def get_ndig(self):
        return self.xaxis.get_ndig()

    @property
    def xdtu_r(self):
        return self.xaxis.coor_r

    @property
    def ydtu_r(self):
        return self.yaxis.coor_r

    @property
    def zdtu_r(self):
        return self.zaxis.coor_r

    @property
    def coor_r(self):
        return [self.xaxis.coor_r, self.yaxis.coor_r, self.zaxis.coor_r]

    @property
    def coor(self):
        return [self.xaxis.coor, self.yaxis.coor, self.zaxis.coor]

    def allclose(self, other, rtol=1e-05, atol=1e-08, equal_nan=False):
        return (self.xaxis.allclose(other.xaxis, rtol=rtol, atol=atol, equal_nan=equal_nan) and
                self.yaxis.allclose(other.yaxis, rtol=rtol, atol=atol, equal_nan=equal_nan) and
                self.zaxis.allclose(other.zaxis, rtol=rtol, atol=atol, equal_nan=equal_nan)
                )

    def closeto(self, other, abserror):
        return self.allclose(other, rtol=0.0, atol=abserror, equal_nan=False)

    @classmethod
    def from_header(cls, hdr):
        """Create a DtuConf object from map-like header"""
        xaxis = DtuAxis.from_header(hdr, name='X')
        yaxis = DtuAxis.from_header(hdr, name='Y')
        zaxis = DtuAxis.from_header(hdr, name='Z')
        dtuconf = DtuConf(xaxis, yaxis, zaxis)
        return dtuconf

    @classmethod
    def from_img(cls, img):
        """Create a DtuConf object from a FITS image"""
        if 'MECS' in img:
            hdr = img['MECS'].header
        else:
            hdr = img[0].header

        return cls.from_header(hdr)

    def __eq__(self, other):
        if isinstance(other, DtuConf):
            return (self.xaxis == other.xaxis and
                    self.yaxis == other.yaxis and
                    self.zaxis == other.zaxis)

        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def describe(self, ndig=None):
        # note: set the number of decimal figures in output to the precision
        # value employed in __eq__() function below
        if ndig is None:
            ndig = self.xaxis.ndig

        output = (
            "<DtuConf instance>\n"
              "- XDTU..: {0.xaxis.coor:{width}.{prec}f}\n"
              "- YDTU..: {0.yaxis.coor:{width}.{prec}f}\n"
              "- ZDTU..: {0.zaxis.coor:{width}.{prec}f}\n"
              "- XDTU_0: {0.xaxis.coor_0:{width}.{prec}f}\n"
              "- YDTU_0: {0.yaxis.coor_0:{width}.{prec}f}\n"
              "- ZDTU_0: {0.zaxis.coor_0:{width}.{prec}f}\n"
              "- XDTU_F: {0.xaxis.coor_f:{width}.{prec}f}\n"
              "- YDTU_F: {0.yaxis.coor_f:{width}.{prec}f}\n"
              "- ZDTU_F: {0.zaxis.coor_f:{width}.{prec}f}\n"
        )
        return output.format(self, width=8, prec=ndig)

    def to_wcs(self):
        """Create a WCS structure for DTU measurements"""
        dtur = self.coor_r
        xfac = dtur[0]
        yfac = -dtur[1]

        dtuwcs = astropy.wcs.WCS(naxis=2)
        dtuwcs.wcs.name = 'DTU WCS'
        dtuwcs.wcs.crpix = [0, 0]
        dtuwcs.wcs.cdelt = [1, 1]
        dtuwcs.wcs.crval = [yfac, xfac]
        dtuwcs.wcs.ctype = ['linear', 'linear']
        dtuwcs.wcs.cunit = ['um', 'um']

        return dtuwcs


def apply_on_axis(func, args):
    """Apply a function to the members of a sequence of DtuAxis"""
    coor_avg = float(func([arg.coor for arg in args]))
    coor_f_avg = float(func([arg.coor_f for arg in args]))
    coor_0_avg = float(func([arg.coor_0 for arg in args]))
    names = [arg.name for arg in args]
    if names:
        name = names[0]
    else:
        # empty list, fallback
        name = "X"
    return DtuAxis(name, coor=coor_avg, coor_f=coor_f_avg, coor_0=coor_0_avg)


def average(*args):
    """Return DtuConf instance with averaged values."""
    xaxis_avg = apply_on_axis(numpy.mean, [arg.xaxis for arg in args])
    yaxis_avg = apply_on_axis(numpy.mean, [arg.yaxis for arg in args])
    zaxis_avg = apply_on_axis(numpy.mean, [arg.zaxis for arg in args])

    return DtuConf(xaxis_avg, yaxis_avg, zaxis_avg)


def maxdiff(*args):
    """Return DtuConf instance with maximum differences."""
    xaxis_avg = apply_on_axis(numpy.ptp, [arg.xaxis for arg in args])
    yaxis_avg = apply_on_axis(numpy.ptp, [arg.yaxis for arg in args])
    zaxis_avg = apply_on_axis(numpy.ptp, [arg.zaxis for arg in args])

    return DtuConf(xaxis_avg, yaxis_avg, zaxis_avg)
