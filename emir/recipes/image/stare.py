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

from numina.recipes import RecipeBase, Parameter, DataProductParameter
from numina.recipes import provides, DataFrame
from numina.recipes.requirements import Requirement

from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import NonLinearityCalibration
from emir.dataproducts import SourcesCatalog

from .shared import DirectImageCommon

_logger = logging.getLogger('numina.recipes.emir')

@provides(DataFrame, SourcesCatalog)
class StareImageRecipe(RecipeBase, DirectImageCommon):
    '''
    The effect of recording images of the sky in a given pointing
    position of the TS


    **Observing modes:**

        * Stare image
    
    '''

    __requires__ = [
        DataProductParameter('master_bpm', MasterBadPixelMask, 
                  'Master bad pixel mask'),       
        DataProductParameter('master_bias', MasterBias, 'Master bias image'),
        DataProductParameter('master_dark', MasterDark, 'Master dark image'),
        DataProductParameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        DataProductParameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        Parameter('extinction', 0.0, 'Mean atmospheric extinction'),
        # FIXME: this parameter is optional 
        Parameter('sources', None, 
                  'List of x, y coordinates to measure FWHM',
                  optional=True),
        Parameter('offsets', None, 'List of pairs of offsets',
                  optional=True),
        Requirement('instrument.detector.channels', 'List of channels'),
        Requirement('instrument.detector.shape', 'Detector shape'),
    ]

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
        
        