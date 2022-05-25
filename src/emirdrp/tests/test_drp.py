
from numina.core import BaseRecipe

from ..loader import load_drp


def test_recipes_are_defined():

    emir_drp = load_drp()

    assert 'default' in emir_drp.pipelines

    for pipeval in emir_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipe = pipeval.get_recipe_object(key)
            assert isinstance(recipe, BaseRecipe)
