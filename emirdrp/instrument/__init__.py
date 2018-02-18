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

"""EMIR constants"""

import math


EMIR_PIXSCALE = 18.0
EMIR_GAIN = 5.0 # ADU / e-
EMIR_RON = 5.69 # ADU
EMIR_NBARS = 55
# Platescale in radians
EMIR_PLATESCALE_RADS = (0.1944 * 1 / 3600.0) * math. pi / 180.0
