#
# Copyright 2013-2014 Universidad Complutense de Madrid
# 
# This file is part of PyEmir
# 
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PyEmir is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
#

'''Tests for AIV Recipes for EMIR'''

import os
from tempfile import mkstemp
import unittest

import pyfits
import numpy as np

from numina.core.oresult import ObservationResult, dataframe_from_list
from numina.core.reciperesult import ErrorRecipeResult
from numina.core.recipeinput import RecipeInput
from numina.core.dataframe import DataFrame

from . import SimpleBiasRecipe

class MMTestCase(unittest.TestCase):

    def test1(self):

        # Create some HDUList in memory
        somefits = []
        nimg = 10
        for i in range(nimg):
            hdu = pyfits.PrimaryHDU(data=np.zeros((10,10), dtype='int16'))
            hdul = pyfits.HDUList([hdu])
            somefits.append(hdul)

         #
        # Create the recipe_input
        obsres = ObservationResult()
        obsres.frames = [dataframe_from_list(val) for val in somefits]
        # empty requirements
        reqs = SimpleBiasRecipe.RecipeRequirements()
        recipe_input = RecipeInput(obsres, reqs)
        recipe = SimpleBiasRecipe()
        #
        result = recipe(recipe_input)
        assert isinstance(result, SimpleBiasRecipe.RecipeResult)

    def test3(self):
        import os
        # Create some HDUList
        somefits = []
        nimg = 10
        for i in range(nimg):
            hdu = pyfits.PrimaryHDU(data=np.zeros((10,10), dtype='int16'))
            hdul = pyfits.HDUList([hdu])
            fd, filename = mkstemp()
            hdul.writeto(filename, clobber=True)
            somefits.append(filename)

        #
        # Create the recipe_input
        obsres = ObservationResult()
        obsres.frames = [dataframe_from_list(val) for val in somefits]
        # empty requirements
        reqs = SimpleBiasRecipe.RecipeRequirements()
        recipe_input = RecipeInput(obsres, reqs)
        #
        recipe = SimpleBiasRecipe()
        result = recipe(recipe_input)

        for fname in somefits:
            os.remove(fname)

        assert isinstance(result, SimpleBiasRecipe.RecipeResult)
