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

from numina.core import BaseRecipe, Parameter, define_requirements, define_result
from numina.core import DataProductRequirement, RecipeRequirements
from numina.core import Product
from numina.logger import log_to_history

from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import TelescopeFocus
from emir.dataproducts import DTUFocus
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import NonLinearityCalibration

__all__ = ['TelescopeRoughFocusRecipe', 
           'TelescopeFineFocusRecipe',
           'DTUFocusRecipe',           
           ]

_logger = logging.getLogger('numina.recipes.emir')

class TelescopeRoughFocusRecipeRequirements(RecipeRequirements):
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask')
    master_bias = DataProductRequirement(MasterBias, 'Master bias image')
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 'Polynomial for non-linearity correction')
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 'Master intensity flatfield')
    objects = Parameter(None, 'List of x-y pair of object coordinates'),
    focus_range = Parameter(None, 'Focus range: begin, end and step')   

class TelescopeRoughFocusRecipeResult(RecipeResult):
    focus = Product(TelescopeFocus)

@define_requirements(TelescopeRoughFocusRecipeRequirements)
@define_result(TelescopeRoughFocusRecipeResult)
class TelescopeRoughFocusRecipe(BaseRecipe):
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

    def __init__(self):
        super(TelescopeRoughFocusRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult, reqs):
        return TelescopeRoughFocusRecipeResult(focus=TelescopeFocus())


class TelescopeFineFocusRecipeRequirements(RecipeRequirements):
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask')
    master_bias = DataProductRequirement(MasterBias, 'Master bias image')
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 'Polynomial for non-linearity correction')
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 'Master intensity flatfield')
    objects = Parameter(None, 'List of x-y pair of object coordinates'),   

class TelescopeFineFocusRecipeResult(RecipeResult):
    focus = Product(TelescopeFocus)

@define_requirements(TelescopeFineFocusRecipeRequirements)
@define_result(TelescopeFineFocusRecipeResult)
class TelescopeFineFocusRecipe(BaseRecipe):
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

    def __init__(self):
        super(TelescopeFineFocusRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult, reqs):
        return TelescopeFineFocusRecipeResult(focus=TelescopeFocus())

class DTUFocusRecipeRequirements(RecipeRequirements):
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask')
    master_bias = DataProductRequirement(MasterBias, 'Master bias image')
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 'Polynomial for non-linearity correction')
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 'Master intensity flatfield')
    objects = Parameter(None, 'List of x-y pair of object coordinates'),   
    msm_pattern = Parameter(None, 'List of x-y pair of slit coordinates'),
    dtu_focus_range = Parameter('dtu_focus_range', None, 'Focus range of the DTU: begin, end and step')
    
class DTUFocusRecipeResult(RecipeResult):
    focus = Product(DTUFocus)

@define_requirements(DTUFocusRecipeRequirements)
@define_result(DTUFocusRecipeResult)
class DTUFocusRecipe(BaseRecipe):
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

    def __init__(self):
        super(DTUFocusRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult, reqs):
        return DTUFocusRecipeResult(focus=DTUFocus())

