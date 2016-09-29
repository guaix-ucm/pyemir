import pytest
import astropy.io.fits as fits
import numpy
import numpy.random

import numina.core
import numina.exceptions

from ..join import JoinDitheredImagesRecipe


def create_wcs():
    import numpy
    from astropy import wcs
    w = wcs.WCS(naxis=2)

    w.wcs.crpix = [50.0, 50.0]
    w.wcs.crval = [0.0, 0.0]
    w.wcs.cdelt = [1.0, 1.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cd = numpy.array([
        [5.41615466e-05, -2.55232165e-07],
        [2.55232165e-07,  5.41615466e-05]
    ])
    return w


def create_frame(val=0, std=100, keys=None):
    size = (100, 100)
    data = numpy.random.normal(val, std, size)
    hdu = fits.PrimaryHDU(data)
    # Update WCS header
    wcs = create_wcs()
    hdu.header = wcs.to_header()
    if keys:
        for k, v in keys.items():
            hdu.header[k] = v

    return fits.HDUList([hdu])


def create_ob(value, nimages, nstare=3, exptime=100.0, starttime=0.0):

    frames = []
    for i in range(nimages):
        base = starttime
        off1 = i * exptime
        off2 = off1 + exptime
        t1 = base + off1
        t2 = base + off2
        keys = {'TSUTC1': t1, 'TSUTC2': t2, 'NUM-NCOM': nstare, 'NUM-SK': 1}
        frame = numina.core.DataFrame(frame=create_frame(val=value, keys=keys))
        frames.append(frame)

    obsresult = numina.core.ObservationResult()
    obsresult.frames = frames
    return obsresult


@pytest.mark.parametrize("nimages,naccum", [(7, 10)])
def test_join(nimages, naccum):

    nstare = 1
    value = 13.0
    starttime = 1030040001.00034
    exptime = 105.0

    obsresult = create_ob(value, nimages, nstare, exptime, starttime)
    obsresult.naccum = 1

    recipe = JoinDitheredImagesRecipe()

    while True:

        rinput = recipe.create_input(
            obresult=obsresult,
        )

        import logging
        #ogging.basicConfig(level=logging.DEBUG)
        result = recipe.run(rinput)
        frame_hdul = result.frame.open()
        assert frame_hdul[0].header['NUM-NCOM'] == nimages * nstare
        accum_hdul = result.accum.open()
        print('frame', obsresult.naccum, frame_hdul[0].data.mean(), frame_hdul[0].data.std(), 2*(1.0/obsresult.naccum))
        print('acuum', obsresult.naccum, accum_hdul[0].data.mean(), accum_hdul[0].data.std(), 2*(1-1.0/obsresult.naccum))
        assert accum_hdul[0].header['NUM-NCOM'] == nimages * nstare * obsresult.naccum

        if obsresult.naccum < naccum:
            # Init next loop
            nobsresult = create_ob(value , nimages, nstare, exptime, starttime)
            nobsresult.naccum = obsresult.naccum + 1
            nobsresult.accum = result.accum
            obsresult = nobsresult
        else:
            break

    frame_hdul = result.frame.open()
    accum_hdul = result.accum.open()

    assert frame_hdul[0].header['NUM-NCOM'] == nimages * nstare
    assert frame_hdul['MAP'].data.max() == nimages
    assert frame_hdul[0].header['TSUTC1'] == starttime
    assert accum_hdul[0].header['NUM-NCOM'] == nimages * nstare * naccum
    assert accum_hdul[0].header['TSUTC1'] == starttime


@pytest.mark.parametrize("nimages", [2, 3, 7])
def test_init_aggregate_result(nimages):
    # Test that works in the first run of the loop
    nstare = 3
    value = 9
    starttime = 1030040001.00034
    exptime = 105.0

    obsresult = create_ob(value, nimages, nstare, exptime, starttime)
    obsresult.naccum = 1
    recipe = JoinDitheredImagesRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
    )

    #import logging
    #logging.basicConfig(level=logging.DEBUG)
    result1 = recipe.run(rinput)

    frame_hdul = result1.frame.open()
    accum_hdul = result1.accum.open()

    assert(len(frame_hdul) == len(accum_hdul))

    for key in ['NUM-NCOM', 'TSUTC1', 'TSUTC2']:
        assert frame_hdul[0].header[key] == accum_hdul[0].header[key]

    assert frame_hdul[0].header['NUM-NCOM'] == nimages * nstare

    assert numpy.allclose(frame_hdul[0].data, accum_hdul[0].data)
