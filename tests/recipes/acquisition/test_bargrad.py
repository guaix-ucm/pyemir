"""Recipe of the mode BAR_POS"""

import pytest
from numina.core import DataFrame
from numina.exceptions import RecipeError

from emirdrp.recipes.acquisition.bargrad import BarDetectionRecipe
from emirdrp.testing.create_ob import create_ob_2


def test_bars_without_keywords(frame_bpm, frame_dark, frame_flat, frame_sky):
    """The images are combined, and the recipe stops without ROTANG and the other keywords"""
    obsresult = create_ob_2(4, 2)
    recipe = BarDetectionRecipe()
    rinput = recipe.create_input(
        obresult=obsresult,
        master_bpm=DataFrame(frame=frame_bpm),
        master_dark=DataFrame(frame=frame_dark),
        master_flat=DataFrame(frame=frame_flat),
        master_sky=DataFrame(frame=frame_sky),
        bars_nominal_positions=[[0, 0, 0, 0, 0]],
    )
    with pytest.raises(RecipeError, match="ROTANG"):
        recipe.run(rinput)
