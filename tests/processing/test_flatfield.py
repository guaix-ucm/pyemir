"""Correction of the flat field"""

import logging

import astropy.io.fits as fits
import numpy

from emirdrp.processing.flatfield import FlatFieldCorrector


def test_flat_nonpositive(caplog):
    """The values <= 0 of the flat are 1, and the flat of the caller is not changed"""
    flat = numpy.full((4, 4), 2.0)
    flat[0, 0] = 0.0
    flat[1, 1] = -3.0
    with caplog.at_level(logging.WARNING, logger="emirdrp.processing.flatfield"):
        corrector = FlatFieldCorrector(flat)
    assert "flat has 2 values <= 0" in caplog.text
    assert flat[0, 0] == 0.0 and flat[1, 1] == -3.0

    img = fits.HDUList([fits.PrimaryHDU(numpy.full((4, 4), 8.0))])
    data = corrector(img)[0].data
    assert data[0, 0] == 8.0 and data[1, 1] == 8.0
    assert data[2, 2] == 4.0


def test_flat_nan(caplog):
    flat = numpy.full((4, 4), 2.0)
    flat[0, 0] = numpy.nan
    with caplog.at_level(logging.WARNING, logger="emirdrp.processing.flatfield"):
        FlatFieldCorrector(flat)
    assert "flat has 1 NaN" in caplog.text
