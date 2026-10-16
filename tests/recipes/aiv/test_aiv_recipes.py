"""Recipes of the AIV modes"""

import numpy
import pytest
from numina.core import DataFrame
from numina.exceptions import RecipeError

from emirdrp.recipes.aiv import TestBiasCorrectRecipe
from emirdrp.recipes.aiv.pointsources import TestPointSourceRecipe
from emirdrp.testing.create_ob import create_ob_2


def test_bias_correct(frame_bpm):
    """The result is the combination of the images, an HDUList"""
    obsresult = create_ob_2(4, 3)
    recipe = TestBiasCorrectRecipe()
    rinput = recipe.create_input(obresult=obsresult, master_bpm=DataFrame(frame=frame_bpm))
    result = recipe.run(rinput)
    hdulist = result.frame.open()
    assert hdulist[0].header["NUM-NCOM"] == 3
    assert numpy.allclose(hdulist[0].data, 4)


def test_point_source_without_keywords(frame_bpm, frame_dark, frame_flat, frame_sky):
    """The images are combined, and the recipe stops without FILTER and the other keywords"""
    obsresult = create_ob_2(4, 2)
    recipe = TestPointSourceRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
        master_bpm=DataFrame(frame=frame_bpm),
        master_dark=DataFrame(frame=frame_dark),
        master_flat=DataFrame(frame=frame_flat),
        master_sky=DataFrame(frame=frame_sky),
    )
    with pytest.raises(RecipeError, match="FILTER"):
        recipe.run(rinput)
