#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

EMIR_BIAS_MODES = ['simple', 'bias', 'single']
EMIR_READ_MODES = ['simple', 'bias', 'single', 'cds', 'fowler', 'ramp']

EMIR_PLATESCALE = 0.1944 # arcsec/pixel

# plate scale, focal plane NASMYTH_A, arcsec / mm
# keyword PSCFP
NASMYTH_A_PS = 1.212 ## arcsec / mm
EMIR_CSU_PS = 1.1746 ## arcsec / mm

# plate scale, detector, arcsec / mm
EMIR_DETECTOR_PS = 10.8
# plate scale, detector, arcsec / pix
EMIR_PS_PIX = EMIR_PLATESCALE
# Pixel size, mm
# PXSZDET
EMIR_PIXSCALE = 18.0 # um
EMIR_PIXSIZE_DETECTOR = EMIR_PIXSCALE * 1e-3

# Pixel projected in FP, mm
# keyword PXSZCSU
EMIR_PIXSIZE_CSU = EMIR_PLATESCALE / EMIR_CSU_PS

# plate scale, detector, arcsec / pix
EMIR_DET_PS_PIX = EMIR_DETECTOR_PS * EMIR_PIXSIZE_DETECTOR

EMIR_GAIN = 5.0 # ADU / e-
EMIR_RON = 5.69 # ADU

EMIR_NBARS = 55
EMIR_NAXIS1 = 2048
EMIR_NAXIS1_ENLARGED = 3400
EMIR_NAXIS2 = 2048

EMIR_NPIXPERSLIT_RECTIFIED = 38
EMIR_MINIMUM_SLITLET_WIDTH_MM = 0.0
EMIR_MAXIMUM_SLITLET_WIDTH_MM = 1000.0
EMIR_VALID_GRISMS = ['J', 'H', 'K', 'LR']
EMIR_VALID_FILTERS = ['J', 'H', 'Ksp', 'YJ', 'HK']

