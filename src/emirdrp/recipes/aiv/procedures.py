#
# Copyright 2014-2025 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""AIV Recipes for EMIR"""

from photutils.geometry.circular_overlap import circular_overlap_grid


def encloses_annulus(x_min, x_max, y_min, y_max, nx, ny, r_in, r_out):
    """Encloses function backported from old photutils"""

    gout = circular_overlap_grid(x_min, x_max, y_min, y_max, nx, ny, r_out, 1, 1)
    gin = circular_overlap_grid(x_min, x_max, y_min, y_max, nx, ny, r_in, 1, 1)
    return gout - gin
