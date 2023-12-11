#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


"""Test the common module"""

import numpy
from numpy.testing import assert_allclose

from ..common import normalize_raw


def test_normalize_raw():

    a = numpy.zeros((10, 10))
    a[:5, :5] = -100
    a[5:, :5] = 10000

    b = normalize_raw(a)

    assert b.min() >= 0
    assert b.max() <= 1


def test_normalize_raw_2():
    tt = numpy.array([-1.0, 655.350, 6553.50, 1e5])
    rtt = numpy.array([0.0, 0.01, 0.1, 1.0])

    ntt = normalize_raw(tt)

    # In this case, min is 0 and max is 1
    assert ntt.min() == 0.0
    assert ntt.max() == 1.0

    assert_allclose(ntt, rtt)
