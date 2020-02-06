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
import six
import numina.exceptions
from numina.core import Requirement, Result, Parameter
import numina.types.array as tarray
from numina.processing.combine import basic_processing_with_combination
from scipy import ndimage
from scipy.ndimage.filters import median_filter
from skimage.feature import canny

import emirdrp.datamodel as datamodel
from emirdrp.core import EMIR_PIXSCALE
from emirdrp.core.recipe import EmirRecipe
import emirdrp.requirements as reqs
import emirdrp.products as prods
from .common import normalize, char_slit
from .common import pinhole_char, pinhole_char2

_logger = logging.getLogger(__name__)


class TestMaskRecipe(EmirRecipe):

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

    median_filter_size = Parameter(5, 'Size of the median box')
    canny_sigma = Parameter(3.0, 'Sigma for the canny algorithm')
    obj_min_size = Parameter(200, 'Minimum size of the slit')
    obj_max_size = Parameter(3000, 'Maximum size of the slit')
    slit_size_ratio = Parameter(4.0, 'Minimum ratio between height and width for slits')

    # Recipe Results
    frame = Result(prods.ProcessedImage)
    positions = Result(tarray.ArrayType)
    positions_alt = Result(tarray.ArrayType)
    slitstable = Result(tarray.ArrayType)
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

        _logger.debug('finding slits')

        # First, prefilter with median
        median_filter_size = rinput.median_filter_size
        canny_sigma = rinput.canny_sigma
        obj_min_size = rinput.obj_min_size
        obj_max_size = rinput.obj_max_size

        data1 = hdulist[0].data
        _logger.debug('Median filter with box %d', median_filter_size)
        data2 = median_filter(data1, size=median_filter_size)

        # Grey level image
        img_grey = normalize(data2)

        # Find edges with canny
        _logger.debug('Find edges with canny, sigma %d', canny_sigma)
        edges = canny(img_grey, sigma=canny_sigma)

        # Fill edges
        _logger.debug('Fill holes')
        fill_slits = ndimage.binary_fill_holes(edges)

        _logger.debug('Label objects')
        label_objects, nb_labels = ndimage.label(fill_slits)
        _logger.debug('%d objects found', nb_labels)
        # Filter on the area of the labeled region
        # Perhaps we could ignore this filtering and
        # do it later?
        _logger.debug('Filter objects by size')
        # Sizes of regions
        sizes = numpy.bincount(label_objects.ravel())

        _logger.debug('Min size is %d', obj_min_size)
        _logger.debug('Max size is %d', obj_max_size)

        mask_sizes = (sizes > obj_min_size) & (sizes < obj_max_size)

        # Filter out regions
        nids, = numpy.where(mask_sizes)

        mm = numpy.in1d(label_objects, nids)
        mm.shape = label_objects.shape

        fill_slits_clean = numpy.where(mm, 1, 0)

        # and relabel
        _logger.debug('Label filtered objects')
        relabel_objects, nb_labels = ndimage.label(fill_slits_clean)
        _logger.debug('%d objects found after filtering', nb_labels)
        ids = list(six.moves.range(1, nb_labels + 1))

        _logger.debug('Find regions and centers')
        regions = ndimage.find_objects(relabel_objects)
        centers = ndimage.center_of_mass(data2, labels=relabel_objects,
                                         index=ids
                                         )

        table = char_slit(data2, regions,
                          slit_size_ratio=rinput.slit_size_ratio
                          )

        result = self.create_result(frame=hdulist,
                                    positions=positions,
                                    positions_alt=positions_alt,
                                    slitstable=table,
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
