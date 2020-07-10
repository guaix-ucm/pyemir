#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Datamodel for EMIR and related functions"""

from __future__ import division

import logging

import numina.datamodel as dm
import astropy.wcs

import emirdrp.instrument


_logger = logging.getLogger(__name__)


class EmirDataModel(dm.DataModel):
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
        'filter': dm.QueryAttribute('filter', str),
        'grism': dm.QueryAttribute('grism', str),
        'insconf': dm.QueryAttribute('insconf', str)
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
            'insconf': ('INSCONF', None),
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

    def do_sky_correction(self, img):
        header = img['primary'].header
        return header.get('SKYADD', True)

    def add_computation_time(self, img, time1, time2):
        img[0].header['NUMUTC1'] = time1.isoformat()
        img[0].header['NUMUTC2'] = time2.isoformat()
        return img


def get_mecs_header(hdulist):
    if 'MECS' in hdulist:
        return hdulist['MECS'].header
    else:
        return hdulist[0].header


def get_dtur_from_img(hdulist):
    mecs_header = get_mecs_header(hdulist)
    return get_dtur_from_header(mecs_header)


def get_dtur_from_header(hdr):
    from emirdrp.instrument.dtuconf import DtuConf
    dtuconf = DtuConf.from_header(hdr)
    return dtuconf.coor, dtuconf.coor_r


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
