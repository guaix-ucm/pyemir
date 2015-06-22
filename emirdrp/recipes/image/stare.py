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

'''
Image mode recipes of EMIR

'''

import logging

from numina.core import Parameter
from numina.core import RecipeRequirements, Product
from numina.core import DataFrameType, define_requirements, define_result
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

_logger = logging.getLogger('numina.recipes.emir')


class StareImageRecipeRequirements(RecipeRequirements):
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


@define_requirements(StareImageRecipeRequirements)
@define_result(StareImageRecipeResult)
class StareImageRecipe(DirectImageCommon):

    '''
    The effect of recording images of the sky in a given pointing
    position of the TS


    **Observing modes:**

        * Stare image

    '''

    def __init__(self):
        super(StareImageRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, ri):

        frame, catalog = self.process(ri,
                                      window=None, subpix=1,
                                      stop_after=DirectImageCommon.PRERED)

        return StareImageRecipeResult(frame=frame, catalog=catalog)
