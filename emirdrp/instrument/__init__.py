#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""EMIR constants"""

import math

# FIXME: duplicated in emirdrp.core

EMIR_PIXSCALE = 18.0
EMIR_GAIN = 5.0 # ADU / e-
EMIR_RON = 5.69 # ADU
EMIR_NBARS = 55
# Platescale in radians

EMIR_PLATESCALE_RADS = (0.1944 * 1 / 3600.0) * math. pi / 180.0
EMIR_REF_IPA = 90.0552 # deg