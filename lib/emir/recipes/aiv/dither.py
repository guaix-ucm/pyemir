#
# Copyright 2014 Universidad Complutense de Madrid
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

'''AIV Recipes for EMIR'''

from __future__ import division

import logging
import math

import numpy
import scipy.interpolate as itpl
import scipy.optimize as opz
from astropy.modeling import models, fitting
from astropy.io import fits
import photutils

from numina.array.recenter import img_box, centering_centroid

from numina import __version__
from numina.core import BaseRecipe, RecipeRequirements, RecipeError
from numina.core import Requirement, Product, DataProductRequirement, Parameter
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement
from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import FlatFieldCorrector, SkyCorrector
from numina.flow import SerialFlow
from numina.flow.processing import DivideByExposure
from numina.flow.node import IdNode
from numina.array import combine

from emir.core import RecipeResult
from emir.core import EMIR_BIAS_MODES
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import DataFrameType, MasterIntensityFlat
from emir.dataproducts import CoordinateList2DType
from emir.dataproducts import ArrayType
from emir.core import gather_info

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"
            
class DitheredImageARecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()

class DitheredImageARecipeResult(RecipeResult):
    frame = Product(DataFrameType)

@define_requirements(DitheredImageARecipeRequirements)
@define_result(DitheredImageARecipeResult)
class DitheredImageARecipeResult(BaseRecipe):

    def __init__(self):
        super(DitheredImageARecipeResult, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
                
        data = numpy.zeros((100, 100))
        hdu = fits.PrimaryHDU(data)

        _logger.debug('update result header')
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdr['IMGOBBL'] = 0
        
        result = DitheredImageARecipeResult(frame=hdulist)
        return result

