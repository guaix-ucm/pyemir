import astropy.io.fits as fits
import numpy
import pytest


@pytest.fixture
def frame_bpm():
    data = numpy.zeros((10, 10), dtype="uint8")
    hdu = fits.PrimaryHDU(data)
    hdu.header["UUID"] = "74263122-4747-4a68-b148-7834d702a733"
    return fits.HDUList([hdu])


@pytest.fixture
def frame_dark():
    data = numpy.zeros((10, 10))
    hdu = fits.PrimaryHDU(data)
    hdu.header["UUID"] = "db4cf732-76a7-4f54-85bf-2db2f0ab46ef"
    return fits.HDUList([hdu])


@pytest.fixture
def frame_flat():
    data = numpy.ones((10, 10))
    hdu = fits.PrimaryHDU(data)
    hdu.header["UUID"] = "f7af05ee-27bd-4a04-957d-28d23b11a532"
    return fits.HDUList([hdu])
