#
# Copyright 2011-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""
Extra functions
"""

import logging
import contextlib

import numpy
import astropy.io.fits as fits

_logger = logging.getLogger(__name__)


def generate_bpm_hdu(hdu):
    """Generate a BPM extension copying an existing HDU"""
    return generate_bpm_data(hdu.data, dtype='uint8', header=hdu.header)


def generate_empty_bpm_hdu(hdu, dtype='uint8'):
    """Generate a BPM extension filled with zeros"""
    data = numpy.zeros_like(hdu.data, dtype)
    return generate_bpm_data(data)


def generate_bpm_data(data, dtype='uint8', header=None):
    """Generate a BPM extension from data."""
    # Mask is assumed to be UINT8
    # casting here
    data = data.astype(dtype)
    hdu_bpm = fits.ImageHDU(
        data,
        header=header
    )

    hdu_bpm.header['EXTNAME'] = 'BPM'
    return hdu_bpm
