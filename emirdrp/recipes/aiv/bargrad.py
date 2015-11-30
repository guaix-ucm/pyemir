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

"""Bar characterization using gradients for EMIR"""

from __future__ import division

import logging

import numpy
from scipy.ndimage.filters import median_filter
from scipy.ndimage import convolve1d
from skimage.filters import threshold_li
import astropy.io.fits as fits

from numina.core import Requirement, Product, Parameter, RecipeError
from numina.core.requirements import ObservationResultRequirement
from numina.core.products import ArrayType
from numina.array.utils import wc_to_pix_1d

from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.products import CoordinateList2DType
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement

from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs
from .common import get_dtur_from_header
from .common import get_cs_from_header, get_csup_from_header

from .bardetect import char_bar_peak_l, char_bar_peak_r

PIXSCALE = 18.0


class BarDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    bars_nominal_positions = Requirement(CoordinateList2DType,
                                         'Nominal positions of the bars'
                                         )
    median_filter_size = Parameter(5, 'Size of the median box')
    average_box_row_size = Parameter(7, 'Number of rows to average for fine centering (odd)')
    average_box_col_size = Parameter(21, 'Number of columns to extract for fine centering (odd)')
    fit_peak_npoints = Parameter(7, 'Number of points to use for fitting the peak (odd)')

    # Recipe Products
    frame = Product(DataFrameType)
    derivative = Product(DataFrameType)
    positions = Product(ArrayType)
    DTU = Product(ArrayType)
    IPA = Product(float)
    csupos = Product(ArrayType)
    csusens = Product(ArrayType)
    param_median_filter_size = Product(int)
    param_average_box_row_size  = Product(int)
    param_average_box_col_size  = Product(int)
    param_fit_peak_npoints  = Product(int)

    def run(self, rinput):

        logger = logging.getLogger('numina.recipes.emir')

        logger.info('starting processing for bars detection')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            ipa = hdr['IPA']
            dtub, dtur = get_dtur_from_header(hdr)
            csupos = get_csup_from_header(hdr)
            csusens = get_cs_from_header(hdr)

        except KeyError as error:
            logger.error(error)
            raise RecipeError(error)

        logger.debug('finding bars')

        arr = hdulist[0].data

        # Median filter
        mfilter_size = rinput.median_filter_size
        logger.debug('median filtering, %d columns', mfilter_size)
        arr_median = median_filter(arr, size=(1, mfilter_size))
        logger.debug('median filtering, %d rows', mfilter_size)
        arr_median = median_filter(arr_median, size=(mfilter_size, 1))

        # Savitsky and Golay (1964) filter to compute the X derivative
        # scipy >= xx has a savgol_filter function
        # for compatibility we do it manually

        # Weights hardcoded for a 9 pixel window

        coeffs_are = [0.066667, 0.05, 0.033333, 0.025, 0, -0.025, -0.03333, -0.0666667]
        arr_deriv = convolve1d(arr_median, coeffs_are, axis=-1)

        positions = []

        xfac = dtur[0] / PIXSCALE
        yfac = -dtur[1] / PIXSCALE

        vec = [yfac, xfac]
        logger.debug('DTU shift is %s', vec)

        # and the table of approx positions of the slits
        barstab = rinput.bars_nominal_positions
        # Currently, we only use fields 0 and 2
        # of the nominal positions file

        # Number or rows used
        # These other parameters cab be tuned also
        bstart = 100
        bend = 1900
        logger.debug('ignoring columns outside %d - %d',bstart, bend-1)

        # Threshold in the derivative image using
        # Li's Minimum Cross Entropy method.
        threshold = threshold_li(numpy.abs(arr_deriv[:,bstart:bend]))
        logger.debug('threshold value in derivative image is %7.2f', threshold)
        # extract a region to average
        wy = (rinput.average_box_row_size // 2)
        wx = (rinput.average_box_col_size // 2)
        logger.debug('extraction window is %d rows, %d cols',2*wy+1, 2*wx+1)
        # Fit the peak with these points
        wfit = 2 * (rinput.fit_peak_npoints // 2) + 1
        logger.debug('fit with %d points', wfit)

        for coords in barstab:
            lbarid = int(coords[0])
            rbarid = lbarid + 55
            ref_y_coor = coords[2] + vec[1]
            prow = wc_to_pix_1d(ref_y_coor) - 1
            fits_row = prow + 1 # FITS pixel index

            logger.debug('looking for bars with ids %d - %d', lbarid, rbarid)
            logger.debug('reference y position is Y %7.2f', ref_y_coor)

            # Find the position of each bar
            xpos, st = char_bar_peak_l(arr_deriv, prow, bstart, bend, threshold, lbarid, wx=wx, wy=wy, wfit=wfit)
            positions.append((lbarid, fits_row, xpos+1, st))

            xpos, st = char_bar_peak_r(arr_deriv, prow, bstart, bend, threshold, rbarid, wx=wx, wy=wy, wfit=wfit)
            positions.append((rbarid, fits_row, xpos+1, st))

        logger.debug('end finding bars')
        result = self.create_result(frame=hdulist,
                                    derivative=fits.PrimaryHDU(data=arr_deriv),
                                    positions=positions,
                                    DTU=dtub,
                                    IPA=ipa,
                                    csupos=csupos,
                                    csusens=csusens,
                                    param_median_filter_size=rinput.median_filter_size,
                                    param_average_box_row_size=2*wy+1,
                                    param_average_box_col_size=2*wx+1,
                                    param_fit_peak_npoints=wfit
                                    )
        return result
