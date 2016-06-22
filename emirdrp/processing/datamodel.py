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
    pass
