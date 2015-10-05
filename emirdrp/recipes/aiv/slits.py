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

"""Recipe to detect slits in the AIV mask"""

from __future__ import division

import logging

import six
import numpy
from scipy import ndimage
from scipy.ndimage.filters import median_filter
from skimage.feature import canny
from numina.core import Product, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.core import RecipeError
from numina.core.products import ArrayType

from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement
from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs
from .common import normalize_raw, char_slit
from .common import get_dtur_from_header

_logger = logging.getLogger('numina.recipes.emir')


class TestSlitDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    median_filter_size = Parameter(5, 'Size of the median box')
    canny_sigma = Parameter(3.0, 'Sigma for the canny algorithm')

    # Recipe Results
    frame = Product(DataFrameType)
    slitstable = Product(ArrayType)
    DTU = Product(ArrayType)
    IPA = Product(float)
    DETPA = Product(float)
    DTUPA = Product(float)

    def run(self, rinput):
        _logger.info('starting slit processing')

        _logger.info('basic image reduction')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            ipa = hdr['IPA']
            detpa = hdr['DETPA']
            dtupa = hdr['DTUPA']
            dtur = get_dtur_from_header(hdr)

        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

        _logger.debug('finding slits')


        # Filter values below 0.0
        _logger.debug('Filter values below 0')
        data1 = hdulist[0].data[:]

        data1[data1 < 0.0] = 0.0
        # First, prefilter with median
        median_filter_size = rinput.median_filter_size
        canny_sigma = rinput.canny_sigma

        _logger.debug('Median filter with box %d', median_filter_size)
        data2 = median_filter(data1, size=median_filter_size)

        # Grey level image
        img_grey = normalize_raw(data2)

        # Find edges with canny
        _logger.debug('Find edges with canny, sigma %d', canny_sigma)
        _logger.debug('Find edges with canny, high threshold %d', canny_sigma)
        edges = canny(img_grey, sigma=canny_sigma,
                      high_threshold=0.1, # approx 6500 counts
                      low_threshold=0.05)
        
        # Fill edges
        _logger.debug('Fill holes')
        # I do a dilation and erosion to fill
        # possible holes in 'edges'
        fill = ndimage.binary_dilation(edges)
        fill2 = ndimage.binary_fill_holes(fill)
        fill_slits = ndimage.binary_erosion(fill2)

        _logger.debug('Label objects')
        label_objects, nb_labels = ndimage.label(fill_slits)
        _logger.debug('%d objects found', nb_labels)
        ids = list(six.moves.range(1, nb_labels + 1))

        _logger.debug('Find regions and centers')
        regions = ndimage.find_objects(label_objects)
        centers = ndimage.center_of_mass(data2, labels=label_objects,
                                         index=ids
                                         )

        table = char_slit(data2, regions,
                          slit_size_ratio=-1.0
                          )

        result = self.create_result(frame=hdulist,
                                    slitstable=table,
                                    DTU=dtur,
                                    IPA=ipa,
                                    DETPA=detpa,
                                    DTUPA=dtupa
                                    )

        return result


class TestSlitMaskDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    median_filter_size = Parameter(4, 'Size of the median box')
    canny_sigma = Parameter(3.0, 'Sigma for the canny algorithm')
    obj_min_size = Parameter(200, 'Minimum size of the slit')
    obj_max_size = Parameter(3000, 'Maximum size of the slit')
    slit_size_ratio = Parameter(4.0, 'Minimum ratio between height and width for slits')

    # Recipe Results
    frame = Product(DataFrameType)
    slitstable = Product(ArrayType)
    DTU = Product(ArrayType)
    IPA = Product(float)
    DETPA = Product(float)
    DTUPA = Product(float)

    def run(self, rinput):
        _logger.info('starting slit processing')

        _logger.info('basic image reduction')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            ipa = hdr['IPA']
            detpa = hdr['DETPA']
            dtupa = hdr['DTUPA']
            dtur = get_dtur_from_header(hdr)

        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

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
        img_grey = normalize_raw(data2)

        # Find edges with canny
        _logger.debug('Find edges with canny, sigma %d', canny_sigma)
        edges = canny(img_grey, sigma=canny_sigma)
        
        # Fill edges
        _logger.debug('Fill holes')
        fill_slits =  ndimage.binary_fill_holes(edges)

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
        #plt.imshow(fill_slits_clean)

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

        result = self.create_result(frame=hdulist, slitstable=table,
                                    DTU=dtur,
                                    IPA=ipa,
                                    DETPA=detpa,
                                    DTUPA=dtupa
                                    )

        return result

