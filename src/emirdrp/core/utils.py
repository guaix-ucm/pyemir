#
# Copyright 2020-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import numpy
import math


def create_rot2d(angle):
    """Create 2D rotation matrix"""
    ca = math.cos(angle)
    sa = math.sin(angle)
    return numpy.array([[ca, -sa], [sa, ca]])
