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

from numina.core import BaseRecipe, Parameter
from numina.core import RecipeInput, RecipeResult, Product
from numina.core import FrameDataProduct, define_input, define_result

from emir.dataproducts import SourcesCatalog

from emir.requirements import MasterBadPixelMask_Requirement, MasterBias_Requirement
from emir.requirements import MasterDark_Requirement, NonLinearityCalibration_Requirement
from emir.requirements import MasterIntensityFlatField_Requirement
from emir.requirements import Channels_Requirement, Extinction_Requirement
from emir.requirements import DetectorShape_Requirement, Offsets_Requirement
from emir.requirements import Catalog_Requirement

from .shared import DirectImageCommon

_logger = logging.getLogger('numina.recipes.emir')

class StareImageRecipeInput(RecipeInput):
    master_bpm = MasterBadPixelMask_Requirement()
    master_bias = MasterBias_Requirement()
    master_dark = MasterDark_Requirement()
    nonlinearity = NonLinearityCalibration_Requirement()
    master_intensity_ff = MasterIntensityFlatField_Requirement()    
    extinction = Extinction_Requirement() 
    sources = Catalog_Requirement()
    offsets = Offsets_Requirement()
    iterations = Parameter(4, 'Iterations of the recipe')
    idc = Channels_Requirement()
    ids = DetectorShape_Requirement()

class StareImageRecipeResult(RecipeResult):
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
        
        frame, catalog = self.process(obresult, baseshape, amplifiers, 
                            offsets=offsets, subpix=1,
                            stop_after=DirectImageCommon.PRERED)
        
        return StareImageRecipeResult(frame=frame, catalog=catalog)
