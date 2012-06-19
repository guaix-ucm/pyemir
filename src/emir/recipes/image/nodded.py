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
Beam switched-nodded image mode recipe of EMIR

'''

from numina.core import BaseRecipe, Parameter, DataProductRequirement
from numina.core import Requirement, RecipeInput, RecipeResult
from numina.core import DataFrame, define_input, define_result
from numina.core import Product, FrameDataProduct 
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import NonLinearityCalibration
from emir.dataproducts import SourcesCatalog

from .shared import DirectImageCommon

class NBImageRecipeInput(RecipeInput):
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask')       
    master_bias = DataProductRequirement(MasterBias, 'Master bias image', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 
              'Polynomial for non-linearity correction')
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 
              'Master intensity flatfield')
    extinction = Parameter(0.0, 'Mean atmospheric extinction') 
    sources = Parameter(None, 'List of x, y coordinates to measure FWHM', optional=True)
    offsets = Parameter(None, 'List of pairs of offsets', optional=True)
    idc = Requirement('List of channels', dest='instrument.detector.channels', hidden=True)
    ids = Requirement('Detector shape', dest='instrument.detector.shape', hidden=True)

class NBImageRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)
    catalog = Product(SourcesCatalog)

@define_input(NBImageRecipeInput)
@define_result(NBImageRecipeResult)
class NBImageRecipe(BaseRecipe, DirectImageCommon):
    '''
    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing the TS in cycles
    between two, or more, sky positions. Displacements are larger
    than the EMIR FOV, so the images have no common area. Used
    for sky subtraction


    **Observing modes:**

        * Nodded/Beamswitched images
    
    '''


    def __init__(self):
        super(NBImageRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        
        baseshape = self.parameters['instrument.detector.shape']
        amplifiers = self.parameters['instrument.detector.channels']
        offsets = self.parameters['offsets']

        frame, catalog = self.process(obresult, baseshape, amplifiers, 
                            offsets=offsets, target_is_sky=False,
                            stop_after=DirectImageCommon.FULLRED)

        return NBImageRecipeResult(frame=frame, catalog=catalog)
