import pytest
import numpy

import numina.core
import numina.exceptions

from emirdrp.recipes.spec.subs import BaseABBARecipe
from emirdrp.recipes.spec.coadd import CoaddABBARecipe
from emirdrp.testing.create_ob import create_ob_internal


@pytest.mark.parametrize("nabba", [1, 2, 3, 4])
def test_coadd(nabba, frame_dark, frame_flat):

    expected_data = numpy.zeros((10, 10))
    expected_data[5] = 4
    expected_data[3] = -4

    results = []
    starttime = 1000000001.00034
    for i in range(nabba):
        pos = [5, 3, 3, 5]
        obsresult_i = create_ob_internal(4, 100, starttime + i * 1000, pos)

        recipe = BaseABBARecipe()
        rinput = recipe.create_input(
            obresult=obsresult_i,
            master_dark=numina.core.DataFrame(frame=frame_dark),
            master_flat=numina.core.DataFrame(frame=frame_flat),
        )
        result = recipe.run(rinput)
        results.append(result)

    obs_result = numina.core.ObservationResult()
    obs_result.frames = [result.reduced_mos_abba for result in results]

    recipe = CoaddABBARecipe()
    rinput = recipe.create_input(
        obresult=obs_result,
    )

    result = recipe.run(rinput)
    hdu = result.reduced_mos_abba.frame[0]
    assert hdu.header["NUM-NCOM"] == nabba * 2
    assert hdu.header["TSUTC1"] == starttime
    assert hdu.header["TSUTC2"] == starttime + (nabba - 1) * 1000 + 4 * 100
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
