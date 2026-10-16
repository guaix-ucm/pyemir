import pytest
import numpy

import numina.exceptions

from emirdrp.recipes.spec.subs import BaseABBARecipe
from emirdrp.testing.create_scenes import create_scene_23
from emirdrp.testing.create_ob import create_ob_abba, exposure_keys, make_frame, make_obresult

STARTTIME = 1000000001.00034


def test_subs():
    size = (100, 100)
    expected_data = numpy.zeros((100, 100))
    expected_data[5] = 4
    expected_data[3] = -4
    # the ABBA sequence, the source in row 5 in A and in row 3 in B
    frames = [
        make_frame(create_scene_23(size, peak=2, pos=pos), keys=exposure_keys(idx, 100.0, STARTTIME))
        for idx, pos in enumerate([5, 3, 3, 5])
    ]
    obs_result = make_obresult(frames)

    recipe = BaseABBARecipe()
    rinput = recipe.create_input(
        obresult=obs_result,
    )
    result = recipe.run(rinput)

    reduced_mos_abba_hdul = result.reduced_mos_abba.frame
    assert reduced_mos_abba_hdul[0].header["TSUTC1"] == STARTTIME
    assert reduced_mos_abba_hdul[0].header["TSUTC2"] == STARTTIME + 400.0

    assert reduced_mos_abba_hdul[0].header["NUM-NCOM"] == 2
    assert numpy.allclose(expected_data, result.reduced_mos_abba.frame[0].data)


def test_subs_raise():
    # only one image, not an ABBA sequence
    frame = make_frame(create_scene_23((100, 100), peak=2, pos=5), keys=exposure_keys(0, 100.0, STARTTIME))
    obsresult = make_obresult([frame])

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

    # import logging
    # logging.basicConfig(level=logging.DEBUG)
    result1 = recipe.run(rinput)

    frame_hdul = result1.reduced_mos_abba.open()
    accum_hdul = result1.accum.open()

    assert len(frame_hdul) == len(accum_hdul)

    for key in ["NUM-NCOM", "TSUTC1", "TSUTC2"]:
        assert frame_hdul[0].header[key] == accum_hdul[0].header[key]

    assert frame_hdul[0].header["NUM-NCOM"] == (nimages / 2) * nstare

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
        assert frame_hdul[0].header["NUM-NCOM"] == nimages * nstare
        accum_hdul = result.accum.open()
        # print('frame', obsresult.naccum, frame_hdul[0].data.mean(),
        # frame_hdul[0].data.std(), 2*(1.0/obsresult.naccum))
        # print('acuum', obsresult.naccum, accum_hdul[0].data.mean(),
        # accum_hdul[0].data.std(), 2*(1-1.0/obsresult.naccum))
        assert accum_hdul[0].header["NUM-NCOM"] == nimages * nstare * obsresult.naccum

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

    assert frame_hdul[0].header["NUM-NCOM"] == nimages * nstare
    assert frame_hdul[0].header["TSUTC1"] == starttime
    assert accum_hdul[0].header["NUM-NCOM"] == nimages * nstare * naccum
    assert accum_hdul[0].header["TSUTC1"] == starttime
