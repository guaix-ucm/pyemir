
import pkg_resources

from numina.core import init_drp_system
from numina.core import import_object
from numina.core import BaseRecipe

import emirdrp.loader as el
from emirdrp.core import EmirRecipe


def test_load_pipeline(monkeypatch):
    """Check that our pipeline can be imported"""

    def mockreturn(group=None):

        ep = pkg_resources.EntryPoint('emir', 'emir.loader')
        monkeypatch.setattr(ep, 'load', lambda: el.load_drp)
        return [ep]

    monkeypatch.setattr(pkg_resources, 'iter_entry_points', mockreturn)

    m = init_drp_system()

    thisdrp = m['EMIR']
    # Import all recipes
    for thispipeline in thisdrp.pipelines.values():
        for key, value in thispipeline.recipes.items():

            recipe = import_object(value)
            assert issubclass(recipe, BaseRecipe)
            # Asume that recipes in emirdrp inherit from EmirRecipe
            if value.startswith('emirdrp'):
                assert issubclass(recipe, EmirRecipe)
