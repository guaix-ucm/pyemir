#
# Copyright 2016-2017 Universidad Complutense de Madrid
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
            'insconf': lambda img: 'EMIR-v1',
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


def get_dtur_from_header(hdr):

    # get DTU things from header
    _logger.info('getting DTU position from header')
    xdtu = hdr['XDTU']
    ydtu = hdr['YDTU']
    zdtu = hdr['ZDTU']

    # Defined even if not in the header
    xdtuf = hdr.get('XDTU_F', 1.0)
    ydtuf = hdr.get('YDTU_F', 1.0)
    xdtu0 = hdr.get('XDTU_0', 0.0)
    ydtu0 = hdr.get('YDTU_0', 0.0)
    _logger.info('XDTU=%6.2f YDTU=%6.2f ZDTU=%6.2f', xdtu, ydtu, zdtu)
    _logger.info('XDTU_F=%6.2f YDTU_F=%6.2f', xdtuf, ydtuf)
    _logger.info('XDTU_0=%6.2f YDTU_0=%6.2f', xdtu0, ydtu0)

    xdtur = (xdtu / xdtuf - xdtu0)
    ydtur = (ydtu / ydtuf - ydtu0)
    _logger.info('XDTU_R=%6.2f YDTU_R=%6.2f', xdtur, ydtur)
    dtu = [xdtu, ydtu, zdtu]
    dtur = [xdtur, ydtur, zdtu]
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
