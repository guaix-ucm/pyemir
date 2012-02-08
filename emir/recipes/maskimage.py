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
Mask image Recipe.

Recipe to process mask images.

**Observing modes:**

  *  Mask imaging

'''

import logging

from numina.recipes import RecipeBase, Parameter, provides

from ..dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from ..dataproducts import MSMPositions
from ..dataproducts import MasterIntensityFlat
from ..dataproducts import NonLinearityCalibration

__all__ = ['Recipe']

_logger = logging.getLogger('emir.recipes')

@provides(MSMPositions)
class Recipe(RecipeBase):
    '''Acquire a target.'''

    __requires__ = [       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        Parameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
    ]

    def __init__(self):
        super(Recipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [MSMPostions()]}