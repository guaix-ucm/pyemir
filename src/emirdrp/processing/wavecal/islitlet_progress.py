#
# Copyright 2018 Universidad Complutense de Madrid
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

"""Wavelength calibration parameters for each grism + filter configuration"""

from __future__ import division, print_function

import sys


def islitlet_progress(islitlet, islitlet_max):
    """Auxiliary function to print out progress in loop of slitlets.

    Parameters
    ----------
    islitlet : int
        Current slitlet number.
    islitlet_max : int
        Maximum slitlet number.

    """
    if islitlet % 10 == 0:
        cout = str(islitlet // 10)
    else:
        cout = '.'
    sys.stdout.write(cout)
    if islitlet == islitlet_max:
        sys.stdout.write('\n')
    sys.stdout.flush()