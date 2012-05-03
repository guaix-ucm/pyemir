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
    Recipes for finding the best focus.

'''

import logging

from numina.recipes import RecipeBase, Parameter, provides
from numina.recipes import DataProductParameter
from numina.logger import log_to_history

from ..dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from ..dataproducts import TelescopeFocus
from ..dataproducts import DTUFocus
from ..dataproducts import MasterIntensityFlat
from ..dataproducts import NonLinearityCalibration

__all__ = ['TelescopeRoughFocusRecipe', 
           'TelescopeFineFocusRecipe',
           'DTUFocusRecipe',           
           ]

_logger = logging.getLogger('emir.recipes')

@provides(TelescopeFocus)
class TelescopeRoughFocusRecipe(RecipeBase):
    '''Recipe to compute the telescope focus.
    
    **Observing modes:**

     * Telescope rough focus
     * Emir focus control 

    **Inputs:**

     * A list of images
     * A list of sky images
     * Bias, dark, flat
     * A model of the detector
     * List of focii 

    **Outputs:**
     * Best focus
    '''
    __requires__ = [       
        DataProductParameter('master_bias', MasterBias, 'Master bias image'),
        DataProductParameter('master_dark', MasterDark, 'Master dark image'),
        DataProductParameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        DataProductParameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        DataProductParameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        Parameter('objects', None, 'List of x-y pair of object coordinates'),
        Parameter('focus_range', None, 'Focus range: begin, end and step')        
    ]

    def __init__(self):
        super(TelescopeRoughFocusRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [TelescopeFocus()]}

@provides(TelescopeFocus)
class TelescopeFineFocusRecipe(RecipeBase):
    '''
    Recipe to compute the telescope focus.
    
    **Observing modes:**

        * Telescope fine focus 

    **Inputs:**

     * A list of images
     * A list of sky images
     * Bias, dark, flat
     * A model of the detector
     * List of focii 

    **Outputs:**
     * Best focus
    
    '''

    __requires__ = [       
        DataProductParameter('master_bias', MasterBias, 'Master bias image'),
        DataProductParameter('master_dark', MasterDark, 'Master dark image'),
        DataProductParameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        DataProductParameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        DataProductParameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        Parameter('objects', None, 'List of x-y pair of object coordinates'),        
    ]

    def __init__(self):
        super(TelescopeFineFocusRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [TelescopeFocus()]}

@provides(DTUFocus)
class DTUFocusRecipe(RecipeBase):
    '''
    Recipe to compute the DTU focus.
    
    **Observing modes:**

        * EMIR focus control

    **Inputs:**

     * A list of images
     * A list of sky images
     * Bias, dark, flat
     * A model of the detector
     * List of focii 

    **Outputs:**
     * Best focus
    
    '''

    __requires__ = [       
        DataProductParameter('master_bias', MasterBias, 'Master bias image'),
        DataProductParameter('master_dark', MasterDark, 'Master dark image'),
        DataProductParameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        DataProductParameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        DataProductParameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        Parameter('objects', None, 'List of x-y pair of object coordinates'),
        Parameter('msm_pattern', None, 'List of x-y pair of slit coordinates'),
        Parameter('dtu_focus_range', None, 
                  'Focus range of the DTU: begin, end and step')        
    ]

    def __init__(self):
        super(DTUFocusRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTUFocus()]}
