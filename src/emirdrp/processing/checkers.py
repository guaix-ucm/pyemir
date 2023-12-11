#
# Copyright 2013-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


from __future__ import division
#
import logging

import numpy
import numina.processing as proc


_logger = logging.getLogger(__name__)


class Checker(proc.Corrector):
    """A Node that checks."""

    def __init__(self):
        super(Checker, self).__init__()

    def run(self, img):
        base = img[0].data
        _logger.debug('running checker after flat')
        # Check NaN and Ceros
        mask1 = ~numpy.isfinite(base)
        if numpy.any(mask1):
            _logger.warning('image has %d NaN', mask1.sum())

        return img
