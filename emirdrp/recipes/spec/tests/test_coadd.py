import pytest
import astropy.io.fits as fits
import numpy

import numina.core
import numina.exceptions

from ..subs import BaseABBARecipe
from ..coadd import CoaddABBARecipe

def create_frame(val=0, pos=3, keys=None):
    data = numpy.zeros((10, 10))
    data[pos] += val
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


@pytest.mark.parametrize("nabba", [1,2,3,4])
def test_coadd(nabba):

    expected_data = numpy.zeros((10, 10))
    expected_data[5] = 4
    expected_data[3] = -4

    results = []
    starttime = 1000000001.00034
    for i in range(nabba):
        pos = [5,3,3,5]
        frames = []
        for j in range(4):
            base = starttime + i * 1000
            off1 = j * 100
            off2 = off1 + 100
            t1 = base + off1
            t2 = base + off2
            keys = {'TSUTC1': t1, 'TSUTC2': t2}
            frame = numina.core.DataFrame(frame=create_frame(val=2, pos=pos[j], keys=keys))
            frames.append(frame)

        obsresult = numina.core.ObservationResult()
        obsresult.frames = frames

        recipe = BaseABBARecipe()
        rinput = recipe.create_input(
            obresult=obsresult,
            master_dark=numina.core.DataFrame(frame=create_frame_dark()),
            master_flat=numina.core.DataFrame(frame=create_frame_flat())
        )
        result = recipe.run(rinput)
        results.append(result)

    obsresult = numina.core.ObservationResult()
    obsresult.frames = [result.reduced_mos_abba for result in results]

    recipe = CoaddABBARecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
    )

    result = recipe.run(rinput)
    hdu = result.reduced_mos_abba.frame[0]
    assert hdu.header['NUM-NCOM'] == nabba * 2
    assert hdu.header['TSUTC1'] == starttime
    assert hdu.header['TSUTC2'] == starttime + (nabba - 1) * 1000 + 4 * 100
    assert numpy.allclose(expected_data, hdu.data)


def test_coadd_empty():

    obsresult = numina.core.ObservationResult()
    obsresult.frames = []

    recipe = CoaddABBARecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
    )

    with pytest.raises(numina.exceptions.RecipeError):
        recipe.run(rinput)
