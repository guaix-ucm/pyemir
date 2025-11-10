import pytest
import numpy

import numina.core
import numina.exceptions

from emirdrp.recipes.image.stare import StareImageBaseRecipe
from emirdrp.testing.create_ob import create_ob_2


@pytest.mark.parametrize("nimages", [1, 2, 3, 4])
def test_stare_with_bpm(nimages, frame_bpm, frame_dark, frame_flat):

    value = 4
    starttime = 1003000001.00034
    exptime = 102.0
    expected_data = value * numpy.ones((10, 10))
    obsresult = create_ob_2(value, nimages, exptime, starttime)

    recipe = StareImageBaseRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
        master_bpm=numina.core.DataFrame(frame=frame_bpm),
        master_dark=numina.core.DataFrame(frame=frame_dark),
        master_flat=numina.core.DataFrame(frame=frame_flat),
    )

    result = recipe.run(rinput)
    hdulist = result.frame.open()

    assert "BPM" in hdulist

    hdu = hdulist[0]
    assert hdu.header["NUM-NCOM"] == nimages
    assert hdu.header["NUM-BPM"] == "uuid:74263122-4747-4a68-b148-7834d702a733"
    assert hdu.header["NUM-DK"] == "uuid:db4cf732-76a7-4f54-85bf-2db2f0ab46ef"
    assert hdu.header["NUM-FF"] == "uuid:f7af05ee-27bd-4a04-957d-28d23b11a532"
    assert hdu.header["TSUTC1"] == starttime
    assert hdu.header["TSUTC2"] == starttime + nimages * exptime
    assert numpy.allclose(expected_data, hdu.data)


@pytest.mark.parametrize("nimages", [1, 2, 3, 4])
def test_stare(nimages, frame_dark, frame_flat):

    value = 9
    starttime = 1030040001.00034
    exptime = 105.0
    expected_data = value * numpy.ones((10, 10))
    obsresult = create_ob_2(value, nimages, exptime, starttime)

    recipe = StareImageBaseRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
        master_dark=numina.core.DataFrame(frame=frame_dark),
        master_flat=numina.core.DataFrame(frame=frame_flat),
    )

    result = recipe.run(rinput)
    hdulist = result.frame.open()

    bpm_ext = 0
    try:
        bpm_ext = hdulist["BPM"]
    except KeyError:
        pass

    assert bpm_ext != 0, "Wrong BPM extension number"

    hdu = hdulist[0]
    assert hdu.header["NUM-NCOM"] == nimages
    assert hdu.header["NUM-DK"] == "uuid:db4cf732-76a7-4f54-85bf-2db2f0ab46ef"
    assert hdu.header["NUM-FF"] == "uuid:f7af05ee-27bd-4a04-957d-28d23b11a532"
    assert hdu.header["TSUTC1"] == starttime
    assert hdu.header["TSUTC2"] == starttime + nimages * exptime
    assert numpy.allclose(expected_data, hdu.data)
