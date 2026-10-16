import numpy

import numina.core

from emirdrp.recipes.spec.stare import StareSpectraRecipe
from emirdrp.testing.create_ob import exposure_keys, make_frame, make_obresult
from emirdrp.testing.create_scenes import create_scene_4847


def test_subs(frame_dark, frame_flat):

    expected_data = 2.5 + numpy.zeros((10, 10))
    starttime = 1000000001.00034
    frames = []
    for idx, (date, value) in enumerate([("2001-10-01T12:33:43.345", 3), ("2001-10-01T12:35:44.345", 2)]):
        keys = exposure_keys(idx, 100.0, starttime)
        keys["DATE-OBS"] = date
        frames.append(make_frame(create_scene_4847(val=value), keys=keys))
    obsresult = make_obresult(frames)

    recipe = StareSpectraRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
        master_dark=numina.core.DataFrame(frame=frame_dark),
        master_flat=numina.core.DataFrame(frame=frame_flat),
    )
    result = recipe.run(rinput)

    stare_hdul = result.stare.open()
    assert stare_hdul[0].header["TSUTC1"] == starttime
    assert stare_hdul[0].header["TSUTC2"] == starttime + 200.0
    assert "NUM-DK" in stare_hdul[0].header
    assert stare_hdul[0].header["NUM-NCOM"] == 2
    assert numpy.allclose(expected_data, stare_hdul[0].data)
