#
# Copyright 2013-2024 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


"""Corrector to remove flat field"""

import datetime
import logging
import sys

import numpy
import numina.processing as proc

if sys.version_info[:2] <= (3, 10):
    datetime.UTC = datetime.timezone.utc


_logger = logging.getLogger(__name__)


class FlatFieldCorrector(proc.Corrector):
    """A Node that corrects a frame from flat-field."""

    def __init__(self, flatdata, datamodel=None, calibid='calibid-unknown', dtype='float32'):

        self.update_variance = False

        super(FlatFieldCorrector, self).__init__(
            datamodel=datamodel,
            calibid=calibid,
            dtype=dtype)

        self.flatdata = flatdata
        self.flatdata[flatdata <= 0] = 1.0  # To avoid NaN
        self.flat_stats = flatdata.mean()

    def run(self, img):
        imgid = self.get_imgid(img)

        _logger.debug('correcting flat in %s', imgid)
        _logger.debug('flat mean is %f', self.flat_stats)

        data = self.datamodel.get_data(img)
        # data = array.correct_flatfield(data, self.flatdata, dtype=self.dtype)
        # Check if flatdata as 0
        mask1 = self.flatdata < 0
        if numpy.any(mask1):
            _logger.warning('flat has %d zeros', mask1.sum())
        mask2 = ~numpy.isfinite(data)
        if numpy.any(mask2):
            _logger.warning('image has %d NaN', mask2.sum())

        result = data / self.flatdata
        result = result.astype(self.dtype)

        # FIXME: not using datamodel
        img['primary'].data = result
        hdr = img['primary'].header
        hdr['NUM-FF'] = self.calibid
        hdr['history'] = 'Flat-field correction with {}'.format(self.calibid)
        hdr['history'] = 'Flat-field correction time {}'.format(
            datetime.datetime.now(datetime.UTC).isoformat())
        hdr['history'] = 'Flat-field correction mean {}'.format(
            self.flat_stats)
        return img
