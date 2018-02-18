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

'''Test the common module'''

import pytest

import numpy
from numpy.testing import assert_allclose

from ..common import normalize_raw


def test_normalize_raw():

    a = numpy.zeros((10, 10))
    a[:5,:5] = -100
    a[5:,:5] = 10000

    b = normalize_raw(a)

    assert b.min() >=0
    assert b.max() <=1

def test_normalize_raw_2():
    tt = numpy.array([-1.0, 655.350, 6553.50, 1e5])
    rtt = numpy.array([0.0, 0.01, 0.1, 1.0])

    ntt = normalize_raw(tt)

    # In this case, min is 0 and max is 1
    assert ntt.min() == 0.0
    assert ntt.max() == 1.0

    assert_allclose(ntt, rtt)

