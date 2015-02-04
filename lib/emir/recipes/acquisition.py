#
# Copyright 2008-2014 Universidad Complutense de Madrid
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
Recipe for the processing of target acquisition images.

**Observing modes:**

    * Target acquisition
'''

import logging

from numina.core import BaseRecipe
from numina.core import RecipeRequirements, Product
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement

from emir.core import RecipeResult
from emir.dataproducts import TelescopeOffset, MSMPositions
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterBadPixelMaskRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement


__all__ = ['TargetAcquisitionRecipe', 'MaskImagingRecipe', 'MaskCheckRecipe']

_logger = logging.getLogger('emir.recipes')


class TargetAcquisitionRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()


class TargetAcquisitionRecipeResult(RecipeResult):
    telescope_offset = Product(TelescopeOffset)


@define_requirements(TargetAcquisitionRecipeRequirements)
@define_result(TargetAcquisitionRecipeResult)
class TargetAcquisitionRecipe(BaseRecipe):

    '''
    Acquire a target.

    Recipe for the processing of target acquisition images.

    **Observing modes:**

        * Target acquisition


    '''

    def __init__(self):
        super(TargetAcquisitionRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult, reqs):
        return self.create_result(telescope_offset=TelescopeOffset())


class MaskImagingRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()


class MaskImagingRecipeResult(RecipeResult):
    msm_positions = Product(MSMPositions)


@define_requirements(MaskImagingRecipeRequirements)
@define_result(MaskImagingRecipeResult)
class MaskImagingRecipe(BaseRecipe):

    '''Acquire a target.

    Mask image Recipe.

    Recipe to process mask images.

    **Observing modes:**

      *  Mask imaging
    '''

    def __init__(self):
        super(MaskImagingRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult, reqs):
        return self.create_result(msm_positions=MSMPositions())


class MaskCheckRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()


class MaskCheckRecipeResult(RecipeResult):
    msm_positions = Product(MSMPositions)
    telescope_offset = Product(TelescopeOffset)


@define_requirements(MaskCheckRecipeRequirements)
@define_result(MaskCheckRecipeResult)
class MaskCheckRecipe(BaseRecipe):

    '''
    Acquire a target.

    Recipe for the processing of multi-slit/long-slit check images.

    **Observing modes:**

        * MSM and LSM check

    '''

    def __init__(self):
        super(MaskCheckRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult, reqs):
        return self.create_result(msm_positions=MSMPositions(),
                                  telescope_offset=TelescopeOffset())
