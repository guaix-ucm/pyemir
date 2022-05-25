
import pytest
import astropy.io.fits as fits
import numpy

import numina.core
import numina.exceptions

from ..subs import BaseABBARecipe

def create_frame_noise(val=0, std=100, pos=3, keys=None):
    size = (100, 100)
    data = numpy.random.normal(val, std, size)
    data[pos] += val
    hdu = fits.PrimaryHDU(data.astype('float32'))
    if keys:
        for k, v in keys.items():
            hdu.header[k] = v

    return fits.HDUList([hdu])


def create_frame(val=0, std=100, pos=3, keys=None):
    size = (100, 100)
    data = numpy.zeros(size)
    data[pos] += val
    hdu = fits.PrimaryHDU(data.astype('float32'))
    if keys:
        for k, v in keys.items():
            hdu.header[k] = v

    return fits.HDUList([hdu])

def create_ob_abba(value, exptime=100.0, starttime=0.0):
    pos = [20, 30, 30, 20]
    frames = []
    nimages = 4
    for i in range(nimages):
        base = starttime
        off1 = i * exptime
        off2 = off1 + exptime
        t1 = base + off1
        t2 = base + off2
        keys = {'TSUTC1': t1, 'TSUTC2': t2, 'NUM-SK': 1}
        frame = numina.core.DataFrame(
            frame=create_frame_noise(
                val=value,
                pos=pos[i],
                keys=keys
            )
        )
        frames.append(frame)

    obsresult = numina.core.ObservationResult()
    obsresult.frames = frames
    obsresult.mode = 'LS_ABBA'
    return obsresult


def test_subs():

    expected_data = numpy.zeros((100, 100))
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

    reduced_mos_abba_hdul = result.reduced_mos_abba.frame
    assert reduced_mos_abba_hdul[0].header['TSUTC1'] == 1000000001.00034
    assert reduced_mos_abba_hdul[0].header['TSUTC2'] == 1000000400.00034

    assert reduced_mos_abba_hdul[0].header['NUM-NCOM'] == 2
    assert numpy.allclose(expected_data, result.reduced_mos_abba.frame[0].data)


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


@pytest.mark.parametrize("nimages", [4])
def test_init_aggregate_result(nimages):
    # Test that works in the first run of the loop
    nstare = 1
    value = 10000
    starttime = 1030040001.00034
    exptime = 105.0

    obsresult = create_ob_abba(value, exptime, starttime)
    obsresult.naccum = 1
    recipe = BaseABBARecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
    )

    #import logging
    #logging.basicConfig(level=logging.DEBUG)
    result1 = recipe.run(rinput)

    frame_hdul = result1.reduced_mos_abba.open()
    accum_hdul = result1.accum.open()

    assert(len(frame_hdul) == len(accum_hdul))

    for key in ['NUM-NCOM', 'TSUTC1', 'TSUTC2']:
        assert frame_hdul[0].header[key] == accum_hdul[0].header[key]

    assert frame_hdul[0].header['NUM-NCOM'] == (nimages / 2) * nstare

    assert numpy.allclose(frame_hdul[0].data, accum_hdul[0].data)


@pytest.mark.parametrize("naccum", [2, 7])
def test_accum_spec(naccum):

    nstare = 1
    nimages = 2
    value = 13.0
    starttime = 1030040001.00034
    exptime = 105.0
    obsresult = create_ob_abba(value, exptime, starttime)
    obsresult.naccum = 1

    recipe = BaseABBARecipe()

    while True:

        rinput = recipe.create_input(
            obresult=obsresult,
        )

        result = recipe.run(rinput)
        frame_hdul = result.reduced_mos_abba.open()
        assert frame_hdul[0].header['NUM-NCOM'] == nimages * nstare
        accum_hdul = result.accum.open()
        # print('frame', obsresult.naccum, frame_hdul[0].data.mean(), frame_hdul[0].data.std(), 2*(1.0/obsresult.naccum))
        # print('acuum', obsresult.naccum, accum_hdul[0].data.mean(), accum_hdul[0].data.std(), 2*(1-1.0/obsresult.naccum))
        assert accum_hdul[0].header['NUM-NCOM'] == nimages * nstare * obsresult.naccum

        if obsresult.naccum < naccum:
            # Init next loop
            nobsresult = create_ob_abba(value, exptime, starttime)
            nobsresult.naccum = obsresult.naccum + 1
            nobsresult.accum = result.accum
            obsresult = nobsresult
        else:
            break

    frame_hdul = result.reduced_mos_abba.open()
    accum_hdul = result.accum.open()

    assert frame_hdul[0].header['NUM-NCOM'] == nimages * nstare
    assert frame_hdul[0].header['TSUTC1'] == starttime
    assert accum_hdul[0].header['NUM-NCOM'] == nimages * nstare * naccum
    assert accum_hdul[0].header['TSUTC1'] == starttime
