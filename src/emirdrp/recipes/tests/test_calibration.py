#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

from emirdrp.loader import load_drp
from emirdrp.recipes.auxiliary import BiasRecipe
import numina.drps

def test_recipe1(drpmocker):

    drpmocker.add_drp('EMIR', load_drp)

    insdrp = numina.drps.get_system_drps().query_by_name('EMIR')
    pipeline = insdrp.pipelines.get('default')
    recipe = pipeline.get_recipe_object('IMAGE_BIAS')

    assert isinstance(recipe, BiasRecipe)
