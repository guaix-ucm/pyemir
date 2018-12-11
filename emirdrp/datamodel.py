#
# Copyright 2016-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyEmir is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
#

"""Datamodel for EMIR and related functions"""

from __future__ import division

import logging

import numina.datamodel
import astropy.wcs

import emirdrp.instrument


_logger = logging.getLogger(__name__)


class QueryAttribute(object):
    def __init__(self, name, tipo, description=""):
        self.name = name
        self.type = tipo
        self.description = description


class EmirDataModel(numina.datamodel.DataModel):
    """Data model of EMIR."""

    meta_dinfo_headers = [
        'readmode',
        'exptime',
        'grism',
        'filter',
        'mode',
        'tstamp',
        'uuid1',
        'uuid2',
        'skyadd'
    ]

    query_attrs = {
        'filter': QueryAttribute('filter', str),
        'grism': QueryAttribute('grism', str),
        'insconf': QueryAttribute('insconf', str)
    }


    def __init__(self):
        defaults = self.default_mappings()
        # FIXME: this should be computed
        defaults['darktime'] = 'exptime'

        instrument_mappings = {
            'readmode': ('READMODE', 'undefined'),
            'grism': ('GRISM', 'undefined'),
            'filter': ('FILTER', 'undefined'),
            'tstamp': ('TSTAMP', 'undefined'),
            'uuid1': ('UUID', 'undefined'),
            'uuid2': ('EMIRUUID', 'undefined'),
            'skyadd': ('SKYADD', True),
            'insconf': lambda img: 'default',
            'blckuuid': lambda img: None,
        }

        defaults.update(instrument_mappings)

        super(EmirDataModel, self).__init__(
            'EMIR',
            defaults
        )

    @property
    def shape(self):
        return (2048, 2048)

    def get_imgid(self, img):
        hdr = self.get_header(img)
        if 'UUID' in hdr:
            return 'uuid:{}'.format(hdr['UUID'])
        elif 'EMIRUUID' in hdr:
            return 'uuid:{}'.format(hdr['EMIRUUID'])
        elif 'TSUTC1' in hdr:
            return 'tsutc:{:16.5f}'.format(hdr['TSUTC1'])
        else:
            return super(EmirDataModel, self).get_imgid(img)

    def do_sky_correction(self, img):
        header = img['primary'].header
        return header.get('SKYADD', True)

    def add_computation_time(self, img, time1, time2):
        img[0].header['NUMUTC1'] = time1.isoformat()
        img[0].header['NUMUTC2'] = time2.isoformat()
        return img

    def get_dtur_from_header(self, hdr):

        return DtuConf.from_header(hdr)


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
        return self.coor / self.coor_f - self.coor_0

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


def get_csup_from_header(hdr):

    # CSUP keys
    default = 0.0
    first_idx = 1
    last_idx = 110
    _logger.info('getting CSUP keys from header')
    values = []
    for idx in range(first_idx, last_idx+1):
        key = "CSUP{}".format(idx)
        values.append(hdr.get(key, default))

    return values


def get_cs_from_header(hdr):

    # CS keys
    default = 0.0
    first_idx = 0
    last_idx = 439
    _logger.info('getting CS keys from header')
    values = []
    for idx in range(first_idx, last_idx+1):
        key = "CS_{}".format(idx)
        values.append(hdr.get(key, default))

    return values


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
