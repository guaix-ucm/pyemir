#
# Copyright 2011-2012 Universidad Complutense de Madrid
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

from numina.core import BaseRecipe, Parameter, DataProductRequirement
from numina.core import DataFrame, RecipeInput, ValidRecipeResult, Product
from numina.core import Requirement, FrameDataProduct, define_input, define_result

from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import NonLinearityCalibration
from emir.dataproducts import SourcesCatalog

from .shared import DirectImageCommon

_logger = logging.getLogger('numina.recipes.emir')

class StareImageRecipeInput(RecipeInput):
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask')       
    master_bias = DataProductRequirement(MasterBias, 'Master bias image', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 
              'Polynomial for non-linearity correction')
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 
              'Master intensity flatfield')
    extinction = Parameter(0.0, 'Mean atmospheric extinction')
    # FIXME: this parameter is optional 
    sources = Parameter(None, 
              'List of x, y coordinates to measure FWHM',
              optional=True)
    offsets = Parameter(None, 'List of pairs of offsets',
              optional=True)
    iterations = Parameter(4, 'Iterations of the recipe')
    idc = Requirement('List of channels', dest='instrument.detector.channels', hidden=True)
    ids = Requirement('Detector shape', dest='instrument.detector.shape', hidden=True)

class StareImageRecipeResult(ValidRecipeResult):
    frame = Product(FrameDataProduct)
    catalog = Product(SourcesCatalog)

@define_input(StareImageRecipeInput)
@define_result(StareImageRecipeResult)    
class StareImageRecipe(BaseRecipe, DirectImageCommon):
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
        
    def run(self, obresult):
                
        baseshape = self.parameters['instrument.detector.shape']
        amplifiers = self.parameters['instrument.detector.channels']
        offsets = self.parameters['offsets']
        
        return self.process(obresult, baseshape, amplifiers, 
                            offsets=offsets, subpix=1,
                            stop_after=DirectImageCommon.PRERED)
        
        
