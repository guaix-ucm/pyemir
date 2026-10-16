import astropy.io.fits as fits
import numpy
import pytest


def calibration(data, uuid):
    """A calibration image with data and the UUID"""
    hdu = fits.PrimaryHDU(data)
    hdu.header["UUID"] = uuid
    return fits.HDUList([hdu])


@pytest.fixture
def frame_bpm():
    return calibration(numpy.zeros((10, 10), dtype="uint8"), "74263122-4747-4a68-b148-7834d702a733")


@pytest.fixture
def frame_dark():
    return calibration(numpy.zeros((10, 10)), "db4cf732-76a7-4f54-85bf-2db2f0ab46ef")


@pytest.fixture
def frame_flat():
    return calibration(numpy.ones((10, 10)), "f7af05ee-27bd-4a04-957d-28d23b11a532")
