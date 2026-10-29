from numina.core import BaseRecipe

from emirdrp.loader import load_drp


def test_recipes_are_defined():

    emir_drp = load_drp()

    assert "default" in emir_drp.pipelines

    for pipeval in emir_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipe = pipeval.get_recipe_object(key)
            assert isinstance(recipe, BaseRecipe)


def test_recipes_log_with_their_module():
    """Each recipe logs with the logger of its module"""
    emir_drp = load_drp()

    for pipeval in emir_drp.pipelines.values():
        for key in pipeval.recipes:
            recipe = pipeval.get_recipe_object(key)
            assert recipe.logger.name == type(recipe).__module__
