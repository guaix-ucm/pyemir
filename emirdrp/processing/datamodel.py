#
# Copyright 2016 Universidad Complutense de Madrid
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

"""Correct Bad Pixels in Image using a BPM image"""

from __future__ import division

import logging
from numina.flow.processing import SimpleDataModel


_logger = logging.getLogger('numina.recipes.emir')


class EmirDataModel(SimpleDataModel):
    """Data model of EMIR.

    Empty for the moment"""

    def __init__(self):
        # Keys
        self._meta = {
            'readmode': ('READMODE', 'undefined'),
            'bunit': ('BUNIT', 'ADU'),
            'texp': ('EXPTIME', None),
            'grism': ('GRISM', 'undefined'),
            'filter': ('FILTER', 'undefined'),
            'obsmode': ('OBSMODE', 'undefined'),
            'tstamp': ('TSTAMP', 'undefined')
        }

    def gather_info_dframe(self, img):
        with img.open() as hdulist:
            info = self.gather_info_hdu(hdulist)
        return info

    def gather_info_hdu(self, hdulist):
        # READMODE is STRING
        meta = {}
        meta['n_ext'] = len(hdulist)
        extnames = [hdu.header.get('extname', '') for hdu in hdulist[1:]]
        meta['name_ext'] = ['PRIMARY'] + extnames
        for key, val in self._meta.items():
            meta[key] = hdulist[0].header.get(val[0], val[1])

        adu_s = False
        if meta['bunit'].lower() == 'adu/s':
            adu_s = True
        meta['adu_s'] = adu_s

        return meta
