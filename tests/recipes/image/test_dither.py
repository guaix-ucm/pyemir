import numpy
import pytest

from emirdrp.recipes.image.join import JoinDitheredImagesRecipe
from emirdrp.testing.create_base import dither_pattern
from emirdrp.testing.create_ob import create_ob_1


@pytest.mark.parametrize("nimages,naccum", [(7, 10)])
def test_join(nimages, naccum):
    nstare = 3
    value = 13.0
    inittime = 1030040001.00034
    starttime = inittime
    exptime = 105.0
    crpix = dither_pattern([50, 50], 0.0, 20.0, nimages)
    obsresult = create_ob_1(value, nimages, crpix, nstare, exptime, starttime)
    obsresult.naccum = 1

    recipe = JoinDitheredImagesRecipe()
    prev_accum = None

    while True:

        rinput = recipe.create_input(obresult=obsresult, accum=prev_accum)

        result = recipe(rinput)
        frame_hdul = result.frame.open()
        assert frame_hdul[0].header["NUM-NCOM"] == nimages * nstare
        accum_hdul = result.accum.open()
        assert accum_hdul[0].header["NUM-NCOM"] == nimages * nstare * obsresult.naccum
        if obsresult.naccum < naccum:
            # Init next loop
            starttime += exptime * (nimages + 1)
            nobsresult = create_ob_1(value, nimages, crpix, nstare, exptime, starttime)
            nobsresult.naccum = obsresult.naccum + 1
            obsresult = nobsresult
            prev_accum = result.accum
        else:
            break

    frame_hdul = result.frame.open()
    accum_hdul = result.accum.open()

    assert frame_hdul[0].header["NUM-NCOM"] == nimages * nstare
    assert frame_hdul["MAP"].data.max() == nimages
    assert frame_hdul[0].header["TSUTC1"] == starttime
    assert frame_hdul[0].header["TSUTC2"] == starttime + nimages * exptime
    assert accum_hdul[0].header["NUM-NCOM"] == nimages * nstare * naccum
    assert accum_hdul[0].header["TSUTC1"] == inittime
    assert accum_hdul[0].header["TSUTC2"] == starttime + nimages * exptime


@pytest.mark.parametrize("nimages", [2, 3, 7])
def test_init_aggregate_result(nimages):
    # Test that works in the first run of the loop
    nstare = 3
    value = 9
    starttime = 1030040001.00034
    exptime = 105.0
    crpix = numpy.array([(50, 50) for _ in range(nimages)])
    obsresult = create_ob_1(value, nimages, crpix, nstare, exptime, starttime)
    obsresult.naccum = 1
    recipe = JoinDitheredImagesRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
    )

    # import logging
    # logging.basicConfig(level=logging.DEBUG)
    result1 = recipe.run(rinput)

    frame_hdul = result1.frame.open()
    accum_hdul = result1.accum.open()

    assert len(frame_hdul) == len(accum_hdul)

    for key in ["NUM-NCOM", "TSUTC1", "TSUTC2"]:
        assert frame_hdul[0].header[key] == accum_hdul[0].header[key]

    assert frame_hdul[0].header["NUM-NCOM"] == nimages * nstare

    assert numpy.allclose(frame_hdul[0].data, accum_hdul[0].data)
