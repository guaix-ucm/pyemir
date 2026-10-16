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


@pytest.mark.parametrize(
    "keys, imgid",
    [
        ({"UUID": "8b3c2d1e-0000-0000-0000-000000000000"}, "uuid:8b3c2d1e-0000-0000-0000-000000000000"),
        ({"EMIRUUID": "8b3c2d1e-0000-0000-0000-000000000000"}, "uuid:8b3c2d1e-0000-0000-0000-000000000000"),
        ({"TSUTC1": 1030040001.00034}, "tsutc:1030040001.00034"),
    ],
)
def test_imgid_with_prefix(keys, imgid):
    """The identifier tells if it is the UUID or the time stamp of the image"""
    hdu = fits.PrimaryHDU()
    for key, value in keys.items():
        hdu.header[key] = value
    assert EmirDataModel().get_imgid(fits.HDUList([hdu])) == imgid
