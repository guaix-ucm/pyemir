import datetime

import astropy.io.fits as fits
import pytest

from emirdrp.datamodel import EmirDataModel


def test_warns_add_computation_time():
    img = fits.HDUList([fits.PrimaryHDU()])
    time1 = datetime.datetime(2026, 10, 13, 12, 0, 0)
    time2 = datetime.datetime(2026, 10, 13, 12, 5, 0)
    with pytest.warns(DeprecationWarning, match="add_computation_time is deprecated"):
        EmirDataModel().add_computation_time(img, time1, time2)
    assert img[0].header["NUMUTC1"] == time1.isoformat()
    assert img[0].header["NUMUTC2"] == time2.isoformat()
