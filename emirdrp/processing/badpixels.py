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
import warnings
import numpy
import scipy.ndimage as snd
from numina.flow.processing import Corrector


_logger = logging.getLogger('numina.recipes.emir')


class BadPixelCorrectorEmir(Corrector):
    """A Node that corrects a frame from bad pixels."""

    def __init__(self, badpixelmask, mark=True, tagger=None,
                 datamodel=None, dtype='float32'):

        warnings.warn("Use numina.flow.processing.BadPixelCorrector",
                      DeprecationWarning)

        super(BadPixelCorrectorEmir, self).__init__(datamodel,
                                                    tagger=tagger,
                                                    dtype=dtype)

        self.bpm = badpixelmask
        self.good_vals = (self.bpm == 0)
        self.bad_vals = ~self.good_vals
        self.median_box = 5

    def _run(self, img):

        # global median
        imgid = self.get_imgid(img)
        _logger.debug('correcting bad pixel mask in %s', imgid)
        base = img[0].data
        interp = base.copy()
        mask1 = ~numpy.isfinite(base)
        # Fill BPM with median values
        median = numpy.median(base[self.good_vals])
        interp[self.bad_vals] = median
        interp = snd.median_filter(interp, size=self.median_box)
        #
        img[0].data[self.bad_vals] = interp[self.bad_vals]

        mask1 = ~numpy.isfinite(base)
        if numpy.any(mask1):
            _logger.warning('image has %d NaN', mask1.sum())

        mask2 = ~numpy.isfinite(img[0].data)
        if numpy.any(mask2):
            _logger.warning('image has %d NaN', mask2.sum())

        return img