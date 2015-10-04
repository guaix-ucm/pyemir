
from numina.core import import_object

from ..loader import load_drp
from numina.core import BaseRecipe


def test_recipes_are_defined():

    thisdrp = load_drp()

    assert thisdrp is not None

    assert 'EMIR' in thisdrp.instruments

    emir_drp = thisdrp.instruments['EMIR']

    assert 'default' in emir_drp.pipelines

    for pipeval in emir_drp.pipelines.values():
        for k, v in pipeval.recipes.items():
            RecipeClass = import_object(v)
            assert issubclass(RecipeClass, BaseRecipe)

