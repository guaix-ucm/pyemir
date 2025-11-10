import pytest
import astropy.io.fits as fits
import numpy

import numina.core
import numina.exceptions

from emirdrp.recipes.spec.stare import StareSpectraRecipe
from emirdrp.testing.create_scenes import create_scene_4847
from emirdrp.testing.create_base import create_image0


def test_subs(frame_dark, frame_flat):

    expected_data = 2.5 + numpy.zeros((10, 10))
    keys = {
        "DATE-OBS": "2001-10-01T12:33:43.345",
        "TSUTC1": 1000000001.00034,
        "TSUTC2": 1000000100.00034,
    }
    scene1 = create_scene_4847(val=3)
    frames1 = numina.core.DataFrame(frame=create_image0(scene1, keys=keys))
    keys = {
        "DATE-OBS": "2001-10-01T12:35:44.345",
        "TSUTC1": 1000000101.00034,
        "TSUTC2": 1000000200.00034,
    }
    scene2 = create_scene_4847(val=2)
    frames2 = numina.core.DataFrame(frame=create_image0(scene2, keys=keys))

    obsresult = numina.core.ObservationResult()
    obsresult.frames = [frames1, frames2]

    recipe = StareSpectraRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
        master_dark=numina.core.DataFrame(frame=frame_dark),
        master_flat=numina.core.DataFrame(frame=frame_flat),
    )
    result = recipe.run(rinput)

    stare_hdul = result.stare.open()
    assert stare_hdul[0].header["TSUTC1"] == 1000000001.00034
    assert stare_hdul[0].header["TSUTC2"] == 1000000200.00034
    assert "NUM-DK" in stare_hdul[0].header
    assert stare_hdul[0].header["NUM-NCOM"] == 2
    assert numpy.allclose(expected_data, stare_hdul[0].data)
