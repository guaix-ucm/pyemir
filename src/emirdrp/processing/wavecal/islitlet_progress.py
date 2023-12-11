#
# Copyright 2018-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
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
