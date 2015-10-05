#
# Copyright 2011-2014 Universidad Complutense de Madrid
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

"""
Image mode recipes of EMIR
"""


from numina.core import Parameter
from numina.core import RecipeInput, Product
from numina.core import DataFrameType, define_input, define_result
from numina.core.requirements import ObservationResultRequirement

from emirdrp.core import RecipeResult
from emirdrp.products import SourcesCatalog

from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import Extinction_Requirement
from emirdrp.requirements import Offsets_Requirement
from emirdrp.requirements import Catalog_Requirement

from .shared import DirectImageCommon


class StareImageRecipeInput(RecipeInput):
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    extinction = Extinction_Requirement()
    sources = Catalog_Requirement()
    offsets = Offsets_Requirement()
    iterations = Parameter(4, 'Iterations of the recipe')


class StareImageRecipeResult(RecipeResult):
    frame = Product(DataFrameType)
    catalog = Product(SourcesCatalog)


@define_input(StareImageRecipeInput)
@define_result(StareImageRecipeResult)
class StareImageRecipe(DirectImageCommon):

    """
    The effect of recording images of the sky in a given pointing
    position of the TS


    **Observing modes:**

        * Stare image

    """

    def run(self, recipe_input):

        frame, catalog = self.process(recipe_input,
                                      window=None, subpix=1,
                                      stop_after=DirectImageCommon.PRERED)

        return StareImageRecipeResult(frame=frame, catalog=catalog)
