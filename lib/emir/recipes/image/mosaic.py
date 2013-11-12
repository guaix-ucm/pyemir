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
Mosaic image mode recipe of EMIR

'''

import logging

from numina.core import BaseRecipe, Parameter, DataProductRequirement
from numina.core import define_requirements, define_result, DataFrame
from numina.core import Requirement, Product, FrameDataProduct, RecipeRequirements

from emir.core import RecipeResult
from emir.dataproducts import SourcesCatalog

_logger = logging.getLogger('numina.recipes.emir')

class MosaicRecipeRequirements(RecipeRequirements):
    # FIXME: this parameter is optional 
    sources = Parameter(None, 'List of x, y coordinates to measure FWHM')

class MosaicRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)
    catalog = Product(SourcesCatalog)

@define_requirements(MosaicRecipeRequirements)
@define_result(MosaicRecipeResult)
class MosaicRecipe(BaseRecipe):
    '''
    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing to a number of
    sky positions, with separations of the order of the EMIR FOV.
    This command is designed to fully cover a given area on the
    sky, but can also be used to point to a number of sky positions
    on which acquisition is only made at the beginning. Supersky
    frame(s) can be built from the image series.

    **Observing modes:**

        * Mosaic images
    
    '''

    def __init__(self):
        super(MosaicRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult, reqs):
        return MosaicRecipeResult(frame=DataFrame(None), catalog=SourcesCatalog())

#
