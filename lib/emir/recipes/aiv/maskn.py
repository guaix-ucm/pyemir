#
# Copyright 2013-2015 Universidad Complutense de Madrid
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
from scipy import ndimage
from scipy.ndimage.filters import median_filter
from skimage.filter import canny

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
from .common import normalize, char_slit
from .common import get_dtur_from_header

_logger = logging.getLogger('numina.recipes.emir')


GAUSS_FWHM_FACTOR = FWHM_G
PIXSCALE = 18.0


class TestMaskRecipe(EmirRecipe):

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

    median_filter_size = Parameter(4, 'Size of the median box')
    canny_sigma = Parameter(3.0, 'Sigma for the canny algorithm')
    obj_min_size = Parameter(200, 'Minimum size of the slit')
    obj_max_size = Parameter(3000, 'Maximum size of the slit')
    slit_size_ratio = Parameter(4.0, 'Minimum ratio between height and width for slits')

    # Recipe Products
    frame = Product(DataFrameType)
    positions = Product(ArrayType)
    positions_alt = Product(ArrayType)
    slitstable = Product(ArrayType)
    DTU = Product(ArrayType)
    filter = Product(str)
    readmode = Product(str)
    IPA = Product(float)
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
        ids = range(1, nb_labels + 1)

        _logger.debug('Find regions and centers')
        regions = ndimage.find_objects(relabel_objects)
        centers = ndimage.center_of_mass(data2, labels=relabel_objects,
                                         index=ids
                                         )

        table = char_slit(data2, regions, centers,
                          slit_size_ratio=rinput.slit_size_ratio
                          )

        result = self.create_result(frame=hdulist,
                                    positions=positions,
                                    positions_alt=positions_alt,
                                    slitstable=table,
                                    filter=filtername,
                                    DTU=[xdtur, ydtur, zdtur],
                                    readmode=readmode,
                                    IPA=ipa,
                                    param_recenter=rinput.recenter,
                                    param_max_recenter_radius=rinput.max_recenter_radius,
                                    param_box_half_size=rinput.box_half_size
                                    )
        return result
