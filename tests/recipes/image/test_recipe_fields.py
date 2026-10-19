"""The requirements and results of the image recipes, declared in the recipes"""

import pytest

from emirdrp.recipes.image import micro, nodded

COMMON = [
    "obresult",
    "master_bpm",
    "master_bias",
    "master_dark",
    "master_flat",
    "extinction",
    "sources",
    "offsets",
]
PHOTOMETRY = ["check_photometry_levels", "check_photometry_actions"]


@pytest.mark.parametrize(
    "recipe, names, alias_in, alias_res",
    [
        (
            nodded.NBImageRecipe,
            COMMON + ["sky_images", "sky_images_sep_time"] + PHOTOMETRY,
            nodded.NBImageRecipeInput,
            nodded.NBImageRecipeResult,
        ),
        (
            micro.MicroditheredImageRecipe,
            COMMON + ["iterations", "sky_images", "sky_images_sep_time"] + PHOTOMETRY + ["subpixelization", "window"],
            micro.MicroditheredImageRecipeInput,
            micro.MicroditheredImageRecipeResult,
        ),
    ],
)
def test_fields(recipe, names, alias_in, alias_res):
    assert list(recipe.RecipeInput.stored()) == names
    assert list(recipe.RecipeResult.stored()) == ["frame", "catalog"]
    # the classes have the names they had as module classes
    assert recipe.RecipeInput.__name__ == f"{recipe.__name__}Input"
    assert alias_in is recipe.RecipeInput
    assert alias_res is recipe.RecipeResult
