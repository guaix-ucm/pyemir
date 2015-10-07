
from numina.core import import_object
from numina.core import BaseRecipe

from ..loader import load_drp
from ..core import EmirRecipe


def test_recipes_are_defined():

    thisdrp = load_drp()

    assert thisdrp is not None

    assert 'EMIR' in thisdrp.instruments

    emir_drp = thisdrp.instruments['EMIR']

    assert 'default' in emir_drp.pipelines

    for pipeval in emir_drp.pipelines.values():
        for key, val in pipeval.recipes.items():
            recipeClass = import_object(val)
            assert issubclass(recipeClass, BaseRecipe)
            # Asume that recipes in emirdrp inherit from EmirRecipe
            if val.startswith('emirdrp'):
                assert issubclass(recipeClass, EmirRecipe)