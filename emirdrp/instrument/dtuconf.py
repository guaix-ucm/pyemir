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

import logging

import astropy.wcs

import emirdrp.instrument


_logger = logging.getLogger(__name__)


class DtuAxis(object):
    def __init__(self, name, coor, coor_f=1.0, coor_0=0.0):
        if name.lower() not in ['x', 'y', 'z']:
            raise ValueError('"name" must be "X", "Y" or "Z')
        self.name = name
        self.coor = coor
        self.coor_f = coor_f
        self.coor_0 = coor_0

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

    def __eq__(self, other):
        # note: set the precision (number of decimal places) to the same
        # number employed in __str__() function above to print out member
        # values
        ndig = 3
        # FIXME: this is eqv to allclose with atol=1e-4?
        if isinstance(other, DtuAxis):
            result = \
                self.name.lower() == other.name.lower() and \
                (round(self.coor, ndig) == round(other.coor, ndig)) and \
                (round(self.coor_f, ndig) == round(other.coor_f, ndig)) and \
                (round(self.coor_0, ndig) == round(other.coor_0, ndig))

            return result
        return NotImplemented

    def __ne__(self, other):
        return not self == other


class DtuConf(object):
    def __init__(self, xaxis, yaxis, zaxis):
        self.xaxis = xaxis
        self.yaxis = yaxis
        self.zaxis = zaxis

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
        xaxis = DtuAxis.from_header(hdr, name='X')
        yaxis = DtuAxis.from_header(hdr, name='Y')
        zaxis = DtuAxis.from_header(hdr, name='Z')
        dtuconf = DtuConf(xaxis, yaxis, zaxis)
        return dtuconf

    def __eq__(self, other):
        if isinstance(other, DtuConf):
            return (self.xaxis == other.xaxis and
                    self.yaxis == other.yaxis and
                    self.zaxis == other.zaxis)

        return NotImplemented

    def __ne__(self, other):
        return not self == other


def get_dtur_from_header(hdr):

    # get DTU things from header
    _logger.info('getting DTU position from header')
    xdtu = hdr['XDTU']
    xdtuf = hdr.get('XDTU_F', 1.0)
    xdtu0 = hdr.get('XDTU_0', 0.0)

    ydtu = hdr['YDTU']
    ydtuf = hdr.get('YDTU_F', 1.0)
    ydtu0 = hdr.get('YDTU_0', 0.0)

    zdtu = hdr['ZDTU']
    zdtuf = hdr.get('ZDTU_F', 1.0)
    zdtu0 = hdr.get('ZDTU_0', 0.0)

    _logger.info('XDTU=%6.2f YDTU=%6.2f ZDTU=%6.2f', xdtu, ydtu, zdtu)
    _logger.info('XDTU_F=%6.2f YDTU_F=%6.2f ZDTU_F=%6.2f', xdtuf, ydtuf, zdtuf)
    _logger.info('XDTU_0=%6.2f YDTU_0=%6.2f ZDTU_F=%6.2f', xdtu0, ydtu0, zdtu0)

    xdtur = (xdtu / xdtuf - xdtu0)
    ydtur = (ydtu / ydtuf - ydtu0)
    zdtur = (zdtu / zdtuf - zdtu0)
    _logger.info('XDTU_R=%6.2f YDTU_R=%6.2f ZDTU_R=%6.2f', xdtur, ydtur, zdtur)
    dtu = [xdtu, ydtu, zdtu]
    dtur = [xdtur, ydtur, zdtur]
    return dtu, dtur


def create_dtu_wcs_header(hdr):

    # get DTU things from header
    xdtu = hdr['XDTU']
    ydtu = hdr['YDTU']

    # Defined even if not in the header
    xdtuf = hdr.get('XDTU_F', 1.0)
    ydtuf = hdr.get('YDTU_F', 1.0)
    xdtu0 = hdr.get('XDTU_0', 0.0)
    ydtu0 = hdr.get('YDTU_0', 0.0)

    xdtur = (xdtu / xdtuf - xdtu0)
    ydtur = (ydtu / ydtuf - ydtu0)

    xfac = xdtur / emirdrp.instrument.EMIR_PIXSCALE
    yfac = -ydtur / emirdrp.instrument.EMIR_PIXSCALE

    # xout = xin + yfac
    # yout = yin + xfac

    dtuwcs = astropy.wcs.WCS(naxis=2)
    dtuwcs.wcs.name = 'DTU WCS'
    dtuwcs.wcs.crpix = [0, 0]
    dtuwcs.wcs.cdelt = [1, 1]
    dtuwcs.wcs.crval = [yfac, xfac]
    dtuwcs.wcs.ctype = ['linear', 'linear']

    return dtuwcs


def create_dtu_wcs_header_um(hdr):

    # get DTU things from header
    xdtu = hdr['XDTU']
    ydtu = hdr['YDTU']

    # Defined even if not in the header
    xdtuf = hdr.get('XDTU_F', 1.0)
    ydtuf = hdr.get('YDTU_F', 1.0)
    xdtu0 = hdr.get('XDTU_0', 0.0)
    ydtu0 = hdr.get('YDTU_0', 0.0)

    xdtur = (xdtu / xdtuf - xdtu0)
    ydtur = (ydtu / ydtuf - ydtu0)

    xfac = xdtur
    yfac = -ydtur

    dtuwcs = astropy.wcs.WCS(naxis=2)
    dtuwcs.wcs.name = 'DTU WCS um'
    dtuwcs.wcs.crpix = [0, 0]
    dtuwcs.wcs.cdelt = [1, 1]
    dtuwcs.wcs.crval = [yfac, xfac]
    dtuwcs.wcs.ctype = ['linear', 'linear']
    dtuwcs.wcs.cunit = ['um', 'um']

    return dtuwcs
