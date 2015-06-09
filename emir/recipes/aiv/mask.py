#
# Copyright 2013-2014 Universidad Complutense de Madrid
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

import numpy

from numina.core import RecipeError
from numina.core import Requirement, Product, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.constants import FWHM_G
from emir.core import EmirRecipe
from emir.dataproducts import DataFrameType
from emir.dataproducts import CoordinateList2DType
from emir.dataproducts import ArrayType
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
from emir.requirements import MasterSkyRequirement

from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs
from .common import pinhole_char, pinhole_char2
from .common import get_dtur_from_header

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"

GAUSS_FWHM_FACTOR = FWHM_G
PIXSCALE = 18.0


class TestPinholeRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    pinhole_nominal_positions = Requirement(CoordinateList2DType,
                                            'Nominal positions of the pinholes'
                                            )
    shift_coordinates = Parameter(True, 'Use header information to'
                                  ' shift the pinhole positions from (0,0) '
                                  'to X_DTU, Y_DTU')
    box_half_size = Parameter(4, 'Half of the computation box size in pixels')
    recenter = Parameter(True, 'Recenter the pinhole coordinates')
    max_recenter_radius = Parameter(2.0, 'Maximum distance for recentering')

    # Recipe Products
    frame = Product(DataFrameType)
    positions = Product(ArrayType)
    positions_alt = Product(ArrayType)
    DTU = Product(ArrayType)
    filter = Product(str)
    readmode = Product(str)
    IPA = Product(float)
    DETPA = Product(float)
    DTUPA = Product(float)
    param_recenter = Product(bool)
    param_max_recenter_radius = Product(float)
    param_box_half_size = Product(float)

    def run(self, rinput):
        _logger.info('starting processing for slit detection')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        _logger.debug('finding pinholes')

        try:
            filtername = hdr['FILTER']
            readmode = hdr['READMODE']
            ipa = hdr['IPA']
            detpa = hdr['DETPA']
            dtupa = hdr['DTUPA']
            xdtur, ydtur, zdtur = get_dtur_from_header(hdr)
        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

        if rinput.shift_coordinates:
            xfac = xdtur / PIXSCALE
            yfac = -ydtur / PIXSCALE

            vec = numpy.array([yfac, xfac])
            _logger.info('shift is %s', vec)
            ncenters = rinput.pinhole_nominal_positions + vec
        else:
            _logger.info('using pinhole coordinates as they are')
            ncenters = rinput.pinhole_nominal_positions

        _logger.info('pinhole characterization')
        positions = pinhole_char(
            hdulist[0].data,
            ncenters,
            box=rinput.box_half_size,
            recenter_pinhole=rinput.recenter,
            maxdist=rinput.max_recenter_radius
        )

        _logger.info('alternate pinhole characterization')
        positions_alt = pinhole_char2(
            hdulist[0].data, ncenters,
            recenter_pinhole=rinput.recenter,
            recenter_half_box=rinput.box_half_size,
            recenter_maxdist=rinput.max_recenter_radius
        )


        result = self.create_result(frame=hdulist,
                                    positions=positions,
                                    positions_alt=positions_alt,
                                    filter=filtername,
                                    DTU=[xdtur, ydtur, zdtur],
                                    readmode=readmode,
                                    IPA=ipa,
                                    DETPA=detpa,
                                    DTUPA=dtupa,
                                    param_recenter=rinput.recenter,
                                    param_max_recenter_radius=rinput.max_recenter_radius,
                                    param_box_half_size=rinput.box_half_size
                                    )
        return result
