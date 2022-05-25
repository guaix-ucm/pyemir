#
# Copyright 2013-2022 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Recipe to detect slits in the AIV mask"""

from __future__ import division


import numpy
import six
from numina.core import Result, Parameter
from numina.exceptions import RecipeError
import numina.types.array as tarray
from numina.processing.combine import basic_processing_with_combination
from scipy import ndimage
from skimage.feature import canny

import emirdrp.datamodel as datamodel
from emirdrp.core.recipe import EmirRecipe
import emirdrp.requirements as reqs
import emirdrp.products as prods
from .common import normalize_raw, char_slit


class TestSlitDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement()

    median_filter_size = Parameter(5, 'Size of the median box')
    canny_sigma = Parameter(3.0, 'Sigma for the canny algorithm')
    canny_high_threshold = Parameter(0.04, 'High threshold for the Canny algorithm')
    canny_low_threshold = Parameter(0.01, 'High threshold for the Canny algorithm')

    # Recipe Results
    frame = Result(prods.ProcessedImage)
    slitstable = Result(tarray.ArrayType)
    DTU = Result(tarray.ArrayType)
    ROTANG = Result(float)
    DETPA = Result(float)
    DTUPA = Result(float)

    def run(self, rinput):
        self.logger.info('starting slit processing')

        self.logger.info('basic image reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            rotang = hdr['ROTANG']
            detpa = hdr['DETPA']
            dtupa = hdr['DTUPA']
            dtub, dtur = datamodel.get_dtur_from_img(hdulist)

        except KeyError as error:
            self.logger.error(error)
            raise RecipeError(error)

        self.logger.debug('finding slits')

        # Filter values below 0.0
        self.logger.debug('Filter values below 0')
        data1 = hdulist[0].data[:]

        data1[data1 < 0.0] = 0.0
        # First, prefilter with median
        median_filter_size = rinput.median_filter_size
        canny_sigma = rinput.canny_sigma

        self.logger.debug('Median filter with box %d', median_filter_size)
        data2 = ndimage.median_filter(data1, size=median_filter_size)

        # Grey level image
        img_grey = normalize_raw(data2)

        # Find edges with Canny
        self.logger.debug('Find edges, Canny sigma %f', canny_sigma)
        # These thresholds corespond roughly with
        # value x (2**16 - 1)
        high_threshold = rinput.canny_high_threshold
        low_threshold = rinput.canny_low_threshold
        self.logger.debug('Find edges, Canny high threshold %f', high_threshold)
        self.logger.debug('Find edges, Canny low threshold %f', low_threshold)
        edges = canny(img_grey, sigma=canny_sigma,
                      high_threshold=high_threshold,
                      low_threshold=low_threshold)
        
        # Fill edges
        self.logger.debug('Fill holes')
        # I do a dilation and erosion to fill
        # possible holes in 'edges'
        fill = ndimage.binary_dilation(edges)
        fill2 = ndimage.binary_fill_holes(fill)
        fill_slits = ndimage.binary_erosion(fill2)

        self.logger.debug('Label objects')
        label_objects, nb_labels = ndimage.label(fill_slits)
        self.logger.debug('%d objects found', nb_labels)
        ids = list(six.moves.range(1, nb_labels + 1))

        self.logger.debug('Find regions and centers')
        regions = ndimage.find_objects(label_objects)
        centers = ndimage.center_of_mass(data2, labels=label_objects,
                                         index=ids
                                         )

        table = char_slit(data2, regions,
                          slit_size_ratio=-1.0
                          )

        result = self.create_result(frame=hdulist,
                                    slitstable=table,
                                    DTU=dtub,
                                    ROTANG=rotang,
                                    DETPA=detpa,
                                    DTUPA=dtupa
                                    )

        return result


class TestSlitMaskDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement()

    median_filter_size = Parameter(5, 'Size of the median box')
    canny_sigma = Parameter(3.0, 'Sigma for the Canny algorithm')
    canny_high_threshold = Parameter(0.04, 'High threshold for the Canny algorithm')
    canny_low_threshold = Parameter(0.01, 'High threshold for the Canny algorithm')
    obj_min_size = Parameter(200, 'Minimum size of the slit')
    obj_max_size = Parameter(3000, 'Maximum size of the slit')
    slit_size_ratio = Parameter(4.0, 'Minimum ratio between height and width for slits')

    # Recipe Results
    frame = Result(prods.DataFrameType)
    slitstable = Result(tarray.ArrayType)
    DTU = Result(tarray.ArrayType)
    ROTANG = Result(float)
    DETPA = Result(float)
    DTUPA = Result(float)

    def run(self, rinput):
        self.logger.info('starting slit processing')

        self.logger.info('basic image reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            rotang = hdr['ROTANG']
            detpa = hdr['DETPA']
            dtupa = hdr['DTUPA']
            dtub, dtur = datamodel.get_dtur_from_img(hdulist)

        except KeyError as error:
            self.logger.error(error)
            raise RecipeError(error)

        self.logger.debug('finding slits')

        # First, prefilter with median
        median_filter_size = rinput.median_filter_size
        canny_sigma = rinput.canny_sigma
        obj_min_size = rinput.obj_min_size
        obj_max_size = rinput.obj_max_size

        data1 = hdulist[0].data
        self.logger.debug('Median filter with box %d', median_filter_size)
        data2 = ndimage.median_filter(data1, size=median_filter_size)

        # Grey level image
        img_grey = normalize_raw(data2)

        # Find edges with Canny
        self.logger.debug('Find edges with Canny, sigma %f', canny_sigma)
        # These thresholds corespond roughly with
        # value x (2**16 - 1)
        high_threshold = rinput.canny_high_threshold
        low_threshold = rinput.canny_low_threshold
        self.logger.debug('Find edges, Canny high threshold %f', high_threshold)
        self.logger.debug('Find edges, Canny low threshold %f', low_threshold)
        edges = canny(img_grey, sigma=canny_sigma,
                      high_threshold=high_threshold,
                      low_threshold=low_threshold)
        # Fill edges
        self.logger.debug('Fill holes')
        fill_slits = ndimage.binary_fill_holes(edges)

        self.logger.debug('Label objects')
        label_objects, nb_labels = ndimage.label(fill_slits)
        self.logger.debug('%d objects found', nb_labels)
        # Filter on the area of the labeled region
        # Perhaps we could ignore this filtering and
        # do it later?
        self.logger.debug('Filter objects by size')
        # Sizes of regions
        sizes = numpy.bincount(label_objects.ravel())

        self.logger.debug('Min size is %d', obj_min_size)
        self.logger.debug('Max size is %d', obj_max_size)

        mask_sizes = (sizes > obj_min_size) & (sizes < obj_max_size)

        # Filter out regions
        nids, = numpy.where(mask_sizes)

        mm = numpy.in1d(label_objects, nids)
        mm.shape = label_objects.shape

        fill_slits_clean = numpy.where(mm, 1, 0)
        # plt.imshow(fill_slits_clean)

        # and relabel
        self.logger.debug('Label filtered objects')
        relabel_objects, nb_labels = ndimage.label(fill_slits_clean)
        self.logger.debug('%d objects found after filtering', nb_labels)
        ids = list(six.moves.range(1, nb_labels + 1))

        self.logger.debug('Find regions and centers')
        regions = ndimage.find_objects(relabel_objects)
        centers = ndimage.center_of_mass(data2, labels=relabel_objects,
                                         index=ids
                                         )

        table = char_slit(data2, regions,
                          slit_size_ratio=rinput.slit_size_ratio
                          )

        result = self.create_result(frame=hdulist, slitstable=table,
                                    DTU=dtub,
                                    ROTANG=rotang,
                                    DETPA=detpa,
                                    DTUPA=dtupa
                                    )

        return result
