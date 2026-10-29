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

import numpy
import numina.processing as proc

_logger = logging.getLogger(__name__)


class FlatFieldCorrector(proc.Corrector):
    """A Node that corrects a frame from flat-field."""

    def __init__(self, flatdata, datamodel=None, calibid="calibid-unknown", dtype="float32"):

        self.update_variance = False

        super(FlatFieldCorrector, self).__init__(datamodel=datamodel, calibid=calibid, dtype=dtype)

        # A copy, with the values <= 0 replaced by 1, to avoid NaN
        nonpositive = flatdata <= 0
        if numpy.any(nonpositive):
            _logger.warning("flat has %d values <= 0, replaced by 1", nonpositive.sum())
        nonfinite = ~numpy.isfinite(flatdata)
        if numpy.any(nonfinite):
            _logger.warning("flat has %d NaN", nonfinite.sum())
        self.flatdata = numpy.where(nonpositive, 1.0, flatdata)
        self.flat_stats = self.flatdata.mean()

    def run(self, img):
        imgid = self.get_imgid(img)

        _logger.debug("correcting flat in %s", imgid)
        _logger.debug("flat mean is %f", self.flat_stats)

        data = self.datamodel.get_data(img)
        # data = array.correct_flatfield(data, self.flatdata, dtype=self.dtype)
        mask2 = ~numpy.isfinite(data)
        if numpy.any(mask2):
            _logger.warning("image has %d NaN", mask2.sum())

        result = data / self.flatdata
        result = result.astype(self.dtype)

        # FIXME: not using datamodel
        img["primary"].data = result
        hdr = img["primary"].header
        hdr["NUM-FF"] = self.calibid
        hdr["history"] = "Flat-field correction with {}".format(self.calibid)
        hdr["history"] = "Flat-field correction time {}".format(datetime.datetime.now(datetime.UTC).isoformat())
        hdr["history"] = "Flat-field correction mean {}".format(self.flat_stats)
        return img
