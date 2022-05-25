
import pytest
import astropy.io.fits as fits
from astropy import wcs
from astropy.modeling.functional_models import Gaussian2D
import numpy
import numpy.random

import numina.core
import numina.exceptions

from ..join import JoinDitheredImagesRecipe


def create_wcs(crpix=None):

    w = wcs.WCS(naxis=2)

    crpix = [50.0, 50.0] if crpix is None else crpix

    w.wcs.crpix = crpix
    w.wcs.crval = [0.0, 0.0]
    w.wcs.cdelt = [1.0, 1.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cd = numpy.array([
        [5.41615466e-05, -2.55232165e-07],
        [2.55232165e-07,  5.41615466e-05]
    ])
    return w


def dither_pattern(center, base_angle, dist, npoints):

    base_angle_rad = base_angle / 180.0 * numpy.pi
    step = 2*numpy.pi / npoints
    angles = base_angle_rad + numpy.arange(0.0, 2 * numpy.pi, step)
    x = center[0] + dist * numpy.cos(angles)
    y = center[1] + dist * numpy.sin(angles)
    return numpy.asarray([x, y]).T


def create_frame(background=0, std=100, crpix=None, keys=None):

    size = (100, 100)
    # Update WCS header
    wcs = create_wcs(crpix)
    y, x = numpy.mgrid[:100, :100]

    data = numpy.random.normal(background, std, size)
    model = Gaussian2D(
        amplitude=30*background,
        x_mean=wcs.wcs.crpix[0]-3.5,
        y_mean=wcs.wcs.crpix[1]-12.8,
        x_stddev=3.0,
        y_stddev=4.0
    )
    data += model(x, y)
    hdu = fits.PrimaryHDU(data)

    hdu.header = wcs.to_header()
    if keys:
        for k, v in keys.items():
            hdu.header[k] = v

    return fits.HDUList([hdu])


def create_ob(value, nimages, crpix, nstare=3, exptime=100.0, starttime=0.0):

    frames = []
    for i in range(nimages):
        base = starttime
        off1 = i * exptime
        off2 = off1 + exptime
        t1 = base + off1
        t2 = base + off2
        keys = {'TSUTC1': t1, 'TSUTC2': t2, 'NUM-NCOM': nstare, 'NUM-SK': 1}
        frame = numina.core.DataFrame(frame=create_frame(background=value, keys=keys, crpix=crpix[i]))
        frames.append(frame)

    obsresult = numina.core.ObservationResult()
    obsresult.frames = frames
    return obsresult


@pytest.mark.parametrize("nimages,naccum", [(7, 10)])
def test_join(nimages, naccum):
    nstare = 3
    value = 13.0
    inittime = 1030040001.00034
    starttime = inittime
    exptime = 105.0
    crpix = dither_pattern([50, 50], 0.0, 20.0, nimages)
    obsresult = create_ob(value, nimages, crpix, nstare, exptime, starttime)
    obsresult.naccum = 1

    recipe = JoinDitheredImagesRecipe()
    prev_accum = None

    while True:

        rinput = recipe.create_input(
            obresult=obsresult,
            accum=prev_accum
        )

        result = recipe(rinput)
        frame_hdul = result.frame.open()
        assert frame_hdul[0].header['NUM-NCOM'] == nimages * nstare
        accum_hdul = result.accum.open()
        assert accum_hdul[0].header['NUM-NCOM'] == nimages * nstare * obsresult.naccum
        if obsresult.naccum < naccum:
            # Init next loop
            starttime += exptime * (nimages + 1)
            nobsresult = create_ob(value , nimages, crpix, nstare, exptime, starttime)
            nobsresult.naccum = obsresult.naccum + 1
            obsresult = nobsresult
            prev_accum = result.accum
        else:
            break

    frame_hdul = result.frame.open()
    accum_hdul = result.accum.open()

    assert frame_hdul[0].header['NUM-NCOM'] == nimages * nstare
    assert frame_hdul['MAP'].data.max() == nimages
    assert frame_hdul[0].header['TSUTC1'] == starttime
    assert frame_hdul[0].header['TSUTC2'] == starttime + nimages * exptime
    assert accum_hdul[0].header['NUM-NCOM'] == nimages * nstare * naccum
    assert accum_hdul[0].header['TSUTC1'] == inittime
    assert accum_hdul[0].header['TSUTC2'] == starttime + nimages * exptime


@pytest.mark.parametrize("nimages", [2, 3, 7])
def test_init_aggregate_result(nimages):
    # Test that works in the first run of the loop
    nstare = 3
    value = 9
    starttime = 1030040001.00034
    exptime = 105.0
    crpix = numpy.array([(50, 50) for _ in range(nimages)])
    obsresult = create_ob(value, nimages, crpix, nstare, exptime, starttime)
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
