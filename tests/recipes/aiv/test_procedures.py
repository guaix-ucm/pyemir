#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


"""Test the procedures module"""

import numpy
from numpy.testing import assert_allclose

from emirdrp.recipes.aiv.procedures import encloses_annulus


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

    aa = encloses_annulus(
        x_min, x_max, y_min, y_max, a.shape[1], a.shape[0], r_in, r_out
    )

    assert_allclose(aa[50, [20, 40, 50, 60]], [0.0, 1.0, 0.0, 1.0])
