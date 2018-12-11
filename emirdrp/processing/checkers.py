#
# Copyright 2013-2018 Universidad Complutense de Madrid
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

from __future__ import division
#
import logging

import numpy
import numina.processing as  proc


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
