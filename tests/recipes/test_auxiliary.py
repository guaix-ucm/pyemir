"""Recipes of the auxiliary modes"""

import pytest
from numina.core import DataFrame
from numina.exceptions import RecipeError

from emirdrp.recipes.auxiliary import DarkRecipe
from emirdrp.testing.create_ob import create_ob_2


def dark_obresult(exptimes):
    obsresult = create_ob_2(4, len(exptimes))
    for frame, exptime in zip(obsresult.frames, exptimes):
        frame.frame[0].header["EXPTIME"] = exptime
        frame.frame[0].header["READMODE"] = "RAMP"
    return obsresult


@pytest.mark.parametrize("exptimes", [[10.0], [10.0, 10.0, 10.0]])
def test_dark(frame_bpm, exptimes):
    recipe = DarkRecipe()
    rinput = recipe.create_input(obresult=dark_obresult(exptimes), master_bpm=DataFrame(frame=frame_bpm))
    hdulist = recipe.run(rinput).darkframe.open()
    assert hdulist[0].header["NUM-NCOM"] == len(exptimes)
    assert hdulist[0].header["CCDMEAN"] == pytest.approx(4)


def test_dark_different_exptime(frame_bpm):
    recipe = DarkRecipe()
    rinput = recipe.create_input(obresult=dark_obresult([10.0, 20.0]), master_bpm=DataFrame(frame=frame_bpm))
    with pytest.raises(RecipeError, match="wrong exposure time"):
        recipe.run(rinput)
