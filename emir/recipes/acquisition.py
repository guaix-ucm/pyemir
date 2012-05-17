#
# Copyright 2008-2012 Universidad Complutense de Madrid
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

from numina.recipes import RecipeBase, Parameter, DataProductParameter
from numina.recipes import provides

from ..dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from ..dataproducts import TelescopeOffset, MSMPositions
from ..dataproducts import MasterIntensityFlat
from ..dataproducts import NonLinearityCalibration

__all__ = ['TargetAcquisitionRecipe', 'MaskImagingRecipe', 'MaskCheckRecipe']

_logger = logging.getLogger('emir.recipes')

@provides(TelescopeOffset)
class TargetAcquisitionRecipe(RecipeBase):
    '''
    Acquire a target.
    
    Recipe for the processing of target acquisition images.

    **Observing modes:**

        * Target acquisition
    
    
    '''

    __requires__ = [       
        DataProductParameter('master_bias', MasterBias, 'Master bias image'),
        DataProductParameter('master_dark', MasterDark, 'Master dark image'),
        DataProductParameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        DataProductParameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        DataProductParameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
    ]

    def __init__(self):
        super(TargetAcquisitionRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [TelescopeOffset()]}
    
@provides(MSMPositions)
class MaskImagingRecipe(RecipeBase):
    '''Acquire a target.
    
    Mask image Recipe.

    Recipe to process mask images.

    **Observing modes:**

      *  Mask imaging
    '''

    __requires__ = [       
        DataProductParameter('master_bias', MasterBias, 'Master bias image'),
        DataProductParameter('master_dark', MasterDark, 'Master dark image'),
        DataProductParameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        DataProductParameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        DataProductParameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
    ]

    def __init__(self):
        super(MaskImagingRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [MSMPositions()]}
    
@provides(TelescopeOffset, MSMPositions)
class MaskCheckRecipe(RecipeBase):
    '''
    Acquire a target.
    
    Recipe for the processing of multi-slit/long-slit check images.

    **Observing modes:**

        * MSM and LSM check 
            
    '''

    __requires__ = [       
        DataProductParameter('master_bias', MasterBias, 'Master bias image'),
        DataProductParameter('master_dark', MasterDark, 'Master dark image'),
        DataProductParameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        DataProductParameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        DataProductParameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
    ]

    def __init__(self):
        super(MaskCheckRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [TelescopeOffset(), MSMPositions()]}