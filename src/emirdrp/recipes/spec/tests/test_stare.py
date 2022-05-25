
import astropy.io.fits as fits
import numpy

import numina.core
import numina.exceptions

from ..stare import StareSpectraRecipe


def create_frame(val=0,  keys=None):
    data = numpy.zeros((10, 10)) + val
    hdu = fits.PrimaryHDU(data)
    if keys:
        for k, v in keys.items():
            hdu.header[k] = v

    return fits.HDUList([hdu])


def create_frame_dark():
    data = numpy.zeros((10, 10))
    hdu = fits.PrimaryHDU(data)
    hdu.header['UUID'] = '11111111111111111111111111111111'
    return fits.HDUList([hdu])


def create_frame_flat():
    data = numpy.ones((10, 10))
    hdu = fits.PrimaryHDU(data)
    hdu.header['UUID'] = '2222222222222222222222222222222'
    return fits.HDUList([hdu])


def test_subs():

    expected_data = 2.5 + numpy.zeros((10, 10))
    keys = {'DATE-OBS': '2001-10-01T12:33:43.345', 'TSUTC1': 1000000001.00034, 'TSUTC2': 1000000100.00034}
    frames1 = numina.core.DataFrame(frame=create_frame(val=3, keys=keys))
    keys = {'DATE-OBS': '2001-10-01T12:35:44.345', 'TSUTC1': 1000000101.00034, 'TSUTC2': 1000000200.00034}
    frames2 = numina.core.DataFrame(frame=create_frame(val=2, keys=keys))

    obsresult = numina.core.ObservationResult()
    obsresult.frames = [frames1, frames2]

    recipe = StareSpectraRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
        master_dark = numina.core.DataFrame(frame=create_frame_dark()),
        master_flat = numina.core.DataFrame(frame=create_frame_flat())
    )
    result = recipe.run(rinput)

    stare_hdul = result.stare.open()
    assert stare_hdul[0].header['TSUTC1'] == 1000000001.00034
    assert stare_hdul[0].header['TSUTC2'] == 1000000200.00034
    assert 'NUM-DK' in stare_hdul[0].header
    assert stare_hdul[0].header['NUM-NCOM'] == 2
    assert numpy.allclose(expected_data, stare_hdul[0].data)
