#
# Copyright 2008-2013 Universidad Complutense de Madrid
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

from numina.core import BaseRecipe, Parameter, DataProductRequirement
from numina.core import RecipeRequirements, Product
from numina.core import define_requirements, define_result

from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import TelescopeOffset, MSMPositions
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import NonLinearityCalibration

__all__ = ['TargetAcquisitionRecipe', 'MaskImagingRecipe', 'MaskCheckRecipe']

_logger = logging.getLogger('emir.recipes')

class TargetAcquisitionRecipeRequirements(RecipeRequirements):
    master_bias = DataProductRequirement(MasterBias, 'Master bias image'),
    master_dark = DataProductRequirement(MasterDark, 'Master dark image'),
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask'),
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 'Master intensity flatfield'),
    

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
        return TargetAcquisitionRecipeResult(telescope_offset=TelescopeOffset())
    
class MaskImagingRecipeRequirements(RecipeRequirements):
    master_bias = DataProductRequirement(MasterBias, 'Master bias image'),
    master_dark = DataProductRequirement(MasterDark, 'Master dark image'),
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask'),
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 'Master intensity flatfield'),
    

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
        return MaskImagingRecipeResult(msm_positions=MSMPositions())
    
class MaskCheckRecipeRequirements(RecipeRequirements):
    master_bias = DataProductRequirement(MasterBias, 'Master bias image'),
    master_dark = DataProductRequirement(MasterDark, 'Master dark image'),
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask'),
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 'Master intensity flatfield'),
    

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
        return MaskCheckRecipeResult(msm_positions=MSMPositions(), 
                                     telescope_offset=TelescopeOffset())
