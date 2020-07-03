#
# Copyright 2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""EMIR constants with units"""

import astropy.units as u

# FIXME: duplicated in emirdrp.core

EMIR_PIXSIZE = 18.0 * u.micron

EMIR_PIXSCALE = 0.1944 * u.arcsec / u.pixel

# Platescale in focal plane of Nasmyth A
GTC_NASMYTH_A_PLATESCALE = 1.212 * u.arcsec / u.mm

# Plate scale, focal plane in the window of the EMIR
# slightly different from NASMYTH_A_PLATESCALE
# keyword PSCFP
EMIR_CSU_PLATESCALE = 1.1746 * u.arcsec / u.mm

# plate scale, detector, arcsec / mm
EMIR_DETECTOR_PLATESCALE = 10.8  * u.arcsec / u.mm

# This value is stored in keyword IPA
EMIR_REF_IPA = 90.0552 * u.deg # deg
