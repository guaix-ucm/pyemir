#
# Copyright 2013-2016 Universidad Complutense de Madrid
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

"""AIV Recipes for EMIR"""

from __future__ import division

import logging

import numpy
from numina.core import RecipeError
from numina.core import Requirement, Product, Parameter
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement

import emirdrp.datamodel as datamodel
from emirdrp.core import EmirRecipe, EMIR_PIXSCALE
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.products import CoordinateList2DType
from emirdrp.products import DataFrameType
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement
from .common import pinhole_char, pinhole_char2

_logger = logging.getLogger('numina.recipes.emir')


class TestPinholeRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
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
    ROTANG = Product(float)
    DETPA = Product(float)
    DTUPA = Product(float)
    param_recenter = Product(bool)
    param_max_recenter_radius = Product(float)
    param_box_half_size = Product(float)

    def run(self, rinput):
        _logger.info('starting processing for slit detection')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        _logger.debug('finding pinholes')

        try:
            filtername = hdr['FILTER']
            readmode = hdr['READMODE']
            rotang = hdr['ROTANG']
            detpa = hdr['DETPA']
            dtupa = hdr['DTUPA']
            dtub, dtur = datamodel.get_dtur_from_header(hdr)
        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

        if rinput.shift_coordinates:
            xdtur, ydtur, zdtur = dtur
            xfac = xdtur / EMIR_PIXSCALE
            yfac = -ydtur / EMIR_PIXSCALE

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
                                    DTU=dtub,
                                    readmode=readmode,
                                    ROTANG=rotang,
                                    DETPA=detpa,
                                    DTUPA=dtupa,
                                    param_recenter=rinput.recenter,
                                    param_max_recenter_radius=rinput.max_recenter_radius,
                                    param_box_half_size=rinput.box_half_size
                                    )
        return result
