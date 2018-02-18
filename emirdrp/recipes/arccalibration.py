#
# Copyright 2011-2016 Universidad Complutense de Madrid
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

"""Auxiliary Recipes for EMIR"""

from __future__ import division, print_function

import logging

import numpy
from numina.core import Requirement, Product, Parameter
from numina.core.products import ArrayType
from numina.core.products import LinesCatalog
from numina.core.requirements import ObservationResultRequirement

from emirdrp.core import EmirRecipe
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.processing.combine import basic_processing_with_combination


_logger = logging.getLogger('numina.recipes.emir')


class ArcCalibrationRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    lines_catalog = Requirement(LinesCatalog, "Catalog of lines")
    polynomial_degree = Parameter(2, 'Polynomial degree of the arc calibration')

    polynomial_coeffs = Product(ArrayType)

    def run(self, rinput):
        _logger.info('starting arc calibration')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        nslits = len(rinput.slits_catalog)
        coeff_table = numpy.zeros((nslits, rinput.polynomial_degree + 1))


        result = self.create_result(polynomial_coeffs=coeff_table)

        return result