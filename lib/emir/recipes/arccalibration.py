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

#from numina.core import BaseRecipeAutoQC
from numina.core import RecipeError
from numina.core import DataFrame#from numina.core import BaseRecipeAutoQC
from numina.core import Requirement, Product, Parameter
from numina.logger import log_to_history
from numina.array.combine import median
# from numina.flow.processing import BadPixelCorrector
from numina.core.requirements import ObservationResultRequirement
from numina.core.requirements import InstrumentConfigurationRequirement

import emir.instrument.channels as allchannels
from emir.core import EMIR_BIAS_MODES
from emir.core import gather_info_frames
from emir.core import EmirRecipe
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import WavelengthCalibration, MasterSpectralFlat
from emir.dataproducts import ChannelLevelStatisticsType
from emir.dataproducts import ChannelLevelStatistics
from emir.dataproducts import SlitTransmissionCalibration
from emir.dataproducts import CoordinateList1DType
from emir.dataproducts import ArrayType
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterBadPixelMaskRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
from emir.requirements import MasterSpectralFlatFieldRequirement
from .aiv.flows import init_filters_bdf
from .aiv.flows import init_filters_bd
from .aiv.flows import init_filters_b
from .aiv.flows import basic_processing_with_combination

_logger = logging.getLogger('numina.recipes.emir')


class ArcCalibrationRecipe(EmirRecipe):
    '''Bla, bla, bla.


    Bla bla bla
    '''

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    arclines_found_positions = Requirement(CoordinateList1DType,
                                           'Nominal positions of the arc lines'
                                          )
    polynomial_degree = Parameter(2, 'Polynomial degree of the arc calibration')

    polynomial_coeffs = Product(ArrayType)

    def run(self, rinput):
        _logger.info('starting arc calibration')

        flow = init_filters_bd(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)


        #
        # RELLENAR
        #

        numpy_array_with_coeff = [0.5, 1.2]
        result = self.create_result(polynomial_coeffs = numpy_array_with_coeff)

        return result
