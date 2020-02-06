#
# Copyright 2013-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""AIV Recipes for EMIR"""

from __future__ import division

import logging

import numpy
import numina.exceptions
from numina.core import Requirement, Result, Parameter
import numina.types.array as tarray
from numina.processing.combine import basic_processing_with_combination

import emirdrp.datamodel as datamodel
from emirdrp.core import EMIR_PIXSCALE
from emirdrp.core.recipe import EmirRecipe
import emirdrp.products as prods
import emirdrp.requirements as reqs
from .common import pinhole_char, pinhole_char2

_logger = logging.getLogger(__name__)


class TestPinholeRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement()

    pinhole_nominal_positions = Requirement(prods.CoordinateList2DType,
                                            'Nominal positions of the pinholes'
                                            )
    shift_coordinates = Parameter(True, 'Use header information to'
                                  ' shift the pinhole positions from (0,0) '
                                  'to X_DTU, Y_DTU')
    box_half_size = Parameter(4, 'Half of the computation box size in pixels')
    recenter = Parameter(True, 'Recenter the pinhole coordinates')
    max_recenter_radius = Parameter(2.0, 'Maximum distance for recentering')

    # Recipe Results
    frame = Result(prods.ProcessedImage)
    positions = Result(tarray.ArrayType)
    positions_alt = Result(tarray.ArrayType)
    DTU = Result(tarray.ArrayType)
    filter = Result(str)
    readmode = Result(str)
    ROTANG = Result(float)
    DETPA = Result(float)
    DTUPA = Result(float)
    param_recenter = Result(bool)
    param_max_recenter_radius = Result(float)
    param_box_half_size = Result(float)

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
            dtub, dtur = datamodel.get_dtur_from_img(hdulist)
        except KeyError as error:
            _logger.error(error)
            raise numina.exceptions.RecipeError(error)

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
