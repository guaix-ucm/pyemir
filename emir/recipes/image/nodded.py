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

from numina.recipes import RecipeBase, Parameter, provides, DataFrame

from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import NonLinearityCalibration
from emir.dataproducts import SourcesCatalog

from .shared import DirectImageCommon
    
@provides(DataFrame, SourcesCatalog)
class NBImageRecipe(RecipeBase, DirectImageCommon):
    '''
    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing the TS in cycles
    between two, or more, sky positions. Displacements are larger
    than the EMIR FOV, so the images have no common area. Used
    for sky subtraction


    **Observing modes:**

        * Nodded/Beamswitched images
    
    '''

    __requires__ = [
        Parameter('master_bpm', MasterBadPixelMask, 
                  'Master bad pixel mask'),       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        Parameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        # FIXME: this parameter is optional 
        Parameter('sources', None, 'List of x, y coordinates to measure FWHM')
    ]

    def __init__(self):
        super(NBImageRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DataFrame(None), SourcesCatalog()]}

