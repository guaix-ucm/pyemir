#
# Copyright 2015 Universidad Complutense de Madrid
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

'''Test the procedures module'''

import pytest

import numpy
from numpy.testing import assert_allclose

from ..procedures import encloses_annulus


def test_encloses_annulus():
    a = numpy.zeros((100, 100))
    r_in = 6
    r_out = 15.5
    xc = 51.4
    yc = 52.9

    x_min = -0.5 - xc
    x_max = a.shape[1] - 0.5 - xc
    y_min = -0.5 - yc
    y_max = a.shape[1] - 0.5 - yc


    aa = encloses_annulus(x_min, x_max, y_min, y_max,
                          a.shape[1], a.shape[0],
                          r_in, r_out)

    assert_allclose(aa[50, [20, 40, 50, 60]], [0.0, 1.0, 0.0, 1.0])
