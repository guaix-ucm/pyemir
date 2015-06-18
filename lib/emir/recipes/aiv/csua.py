#
# Copyright 2015 Universidad Complutense de Madrid
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

'''Recipe to check the alignment of the CSU mask'''

from __future__ import division

import logging

import numpy
from scipy import ndimage
from scipy.ndimage.filters import median_filter
from skimage.filter import canny

from numina.core import Product, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.core import RecipeError
#
from emir.core import EmirRecipe
from emir.dataproducts import DataFrameType
from emir.dataproducts import ArrayType
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
from emir.requirements import MasterSkyRequirement

from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs
from .common import normalize, char_slit
from .common import get_dtur_from_header

_logger = logging.getLogger('numina.recipes.emir')


class TestCSUAlignmentRecipe(EmirRecipe):

    # Recipe Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    median_filter_size = Parameter(4, 'Size of the median box')

    # Recipe Results
    frame = Product(DataFrameType)
    slitstable = Product(ArrayType)

    def run(self, rinput):
        _logger.info('starting CSU alignment processing')

        _logger.info('basic image reduction')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        _logger.debug('finding slits')

        # First, prefilter with median
        median_filter_size = rinput.median_filter_size
        
        
        data1 = hdulist[0].data
        _logger.debug('Median filter with box %d', median_filter_size)
        data2 = median_filter(data1, size=median_filter_size)

        result = self.create_result(frame=hdulist, slitstable=[1])

        return result

