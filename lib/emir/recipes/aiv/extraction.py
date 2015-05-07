#
# Copyright 2011-2015 Universidad Complutense de Madrid
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

'''Auxiliary Recipes for EMIR'''

import logging

import numpy
from scipy.ndimage.filters import median_filter

#from numina.core import BaseRecipeAutoQC
from numina.core import RecipeError
from numina.core import DataFrame#from numina.core import BaseRecipeAutoQC
from numina.core import Requirement, Product, Parameter
from numina.logger import log_to_history
from numina.array.combine import median
# from numina.flow.processing import BadPixelCorrector
from numina.core.requirements import ObservationResultRequirement

from emir.core import gather_info_frames
from emir.core import EmirRecipe
from emir.dataproducts import SlitTransmissionCalibration
from emir.dataproducts import CoordinateList1DType
from emir.dataproducts import ArrayType
from emir.dataproducts import DataFrameType

from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
from emir.requirements import MasterSkyRequirement
from .flows import init_filters_bdfs
from .flows import basic_processing_with_combination

_logger = logging.getLogger('numina.recipes.emir')


class MaskSpectraExtractionRecipe(EmirRecipe):
    '''
    '''

    # Recipe Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    median_filter_size = Parameter(5, 'Size of the median box')

    frame = Product(DataFrameType)
    #slitstable = Product(ArrayType)
    #DTU = Product(ArrayType)
    #IPA = Product(float)
    #DETPA = Product(float)
    #DTUPA = Product(float)


    def run(self, rinput):
        _logger.info('starting arc calibration')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        # First, prefilter with median
        median_filter_size = rinput.median_filter_size
        
        
        data1 = hdulist[0].data
        _logger.debug('Median filter with box %d', median_filter_size)
        data2 = median_filter(data1, size=median_filter_size)

        # Extract things here

        #
        # RELLENAR
        #
        result = self.create_result(frame=hdulist)

        return result
