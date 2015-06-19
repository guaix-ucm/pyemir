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

from numina.core import Product, Parameter, Requirement
from numina.core.requirements import ObservationResultRequirement
from numina.core import RecipeError
from numina.array.utils import image_box2d
#
from emir.core import EmirRecipe
from emir.dataproducts import CoordinateList2DType
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


def slit_alignment(img, x, y, canny_sigma=5):
    
    # size of the cut box used to analyze the slit
    hbox_x = 50
    hbox_y = 75
    
    # Remove nstrip rows on top and bottom
    # before computing the center of mass
    nstrip = 30
    
    # Fill bfill columns before doing the binary fill
    bfill = 2
    
    # Region of a box, 50 *2  pixels in X, 75 *2 pixels in Y
    roi = image_box2d(x, y, img.shape, (hbox_y, hbox_x))
    
    mm = img[roi]
    
    edges = canny(mm, sigma=canny_sigma)
    
    # Fill edges
    # Fill columns in the borders to help closing the bars
    edges[:,:bfill] = 1
    edges[:,-bfill:] = 1
    fill_slits =  ndimage.binary_fill_holes(edges)
    
    # Reverse image
    vv = (mm) * (1 - fill_slits)

    # Cut nstrip rows
    vv2 = vv[nstrip:-nstrip]
    
    # Center of mass in X
    XX = numpy.arange(vv2.shape[1])
    YY = numpy.arange(vv2.shape[0])
    
    SF1 = vv2.sum(axis=1)
    XB = (vv2 * XX).sum(axis=1) / SF1
    
    SF = SF1.sum()
    XM = (vv2 * XX).sum() / SF
    YM = (vv2 * YY[:,None]).sum() / SF # Broadcasting here YY

    bb = numpy.polyfit(YY, XB, deg=1)
    
    XM += roi[1].start
    YM += roi[0].start + nstrip
    
    return [XM, YM], bb


class TestCSUAlignmentRecipe(EmirRecipe):

    # Recipe Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    slits_positions = Requirement(CoordinateList2DType,
                                  'Positions of the center of the slits'
                                 )

    median_filter_size = Parameter(5, 'Size of the median box')
    canny_sigma = Parameter(5.0, 'Sigma for the canny algorithm')

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

        # First, prefilter with median
        median_filter_size = rinput.median_filter_size
        
        data1 = hdulist[0].data
        _logger.debug('Median filter with box %d', median_filter_size)
        data2 = median_filter(data1, size=median_filter_size)

        nresults = 4
        slits_t = numpy.empty((rinput.slits_positions.shape[0], nresults))
        # Work with slits
        for idx, (x, y) in enumerate(rinput.slits_positions):
            _logger.debug('processing slit in %f %f', x-1, y-1)
            centroid, fit = slit_alignment(data2, x - 1, y -1,
                                           rinput.canny_sigma)
            slits_t[idx, :2] = centroid
            slits_t[idx, 2:] = fit

        # FITS coordinates
        slits_t[:,:2] += 1

        result = self.create_result(frame=hdulist, slitstable=slits_t)

        return result

