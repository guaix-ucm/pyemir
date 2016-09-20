
import pytest
import astropy.io.fits as fits
import numpy

import numina.core
import numina.exceptions

from ..subs import BaseABBARecipe

def create_frame(val=0, pos=3, keys=None):
    data = numpy.zeros((10, 10))
    data[pos] += val
    hdu = fits.PrimaryHDU(data)
    if keys:
        for k, v in keys.items():
            hdu.header[k] = v

    return fits.HDUList([hdu])


def test_subs():

    expected_data = numpy.zeros((10, 10))
    expected_data[5] = 4
    expected_data[3] = -4
    keys = {'TSUTC1': 1000000001.00034, 'TSUTC2': 1000000100.00034}
    frames1 = numina.core.DataFrame(frame=create_frame(val=2, pos=5, keys=keys))
    keys = {'TSUTC1': 1000000101.00034, 'TSUTC2': 1000000200.00034}
    frames2 = numina.core.DataFrame(frame=create_frame(val=2, pos=3, keys=keys))
    keys = {'TSUTC1': 1000000201.00034, 'TSUTC2': 1000000300.00034}
    frames3 = numina.core.DataFrame(frame=create_frame(val=2, pos=3, keys=keys))
    keys = {'TSUTC1': 1000000301.00034, 'TSUTC2': 1000000400.00034}
    frames4 = numina.core.DataFrame(frame=create_frame(val=2, pos=5, keys=keys))

    obsresult = numina.core.ObservationResult()
    obsresult.frames = [frames1, frames2, frames3, frames4]

    recipe = BaseABBARecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
    )
    result = recipe.run(rinput)

    spec_abba_hdul = result.spec_abba.frame
    assert spec_abba_hdul[0].header['TSUTC1'] == 1000000001.00034
    assert spec_abba_hdul[0].header['TSUTC2'] == 1000000400.00034

    assert spec_abba_hdul[0].header['NUM-NCOM'] == 2
    assert numpy.allclose(expected_data, result.spec_abba.frame[0].data)


def test_subs_raise():

    expected_data = numpy.zeros((10, 10))
    expected_data[5] = 4
    expected_data[3] = -4
    keys = {'TSUTC1': 1000000001.00034, 'TSUTC2': 1000000100.00034}
    frames1 = numina.core.DataFrame(frame=create_frame(val=2, pos=5, keys=keys))

    obsresult = numina.core.ObservationResult()
    obsresult.frames = [frames1]

    recipe = BaseABBARecipe()
    rinput = BaseABBARecipe.RecipeInput(
        obresult=obsresult,
    )

    with pytest.raises(numina.exceptions.RecipeError):
        recipe.run(rinput)
