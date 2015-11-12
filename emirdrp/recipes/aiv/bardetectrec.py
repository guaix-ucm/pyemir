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

"""Bar detection recipe for EMIR"""

from __future__ import division

import logging

from scipy.ndimage.filters import median_filter
from skimage.feature import canny

from numina.core import Requirement, Product, Parameter, RecipeError
from numina.core.requirements import ObservationResultRequirement
from numina.core.products import ArrayType

from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.products import CoordinateList2DType
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement

from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs
from .common import normalize_raw
from .common import get_dtur_from_header
from .common import create_dtu_wcs_header
from .common import create_dtu_wcs_header_um
from .bardetect import find_position
from .bardetect import locate_bar_l, locate_bar_r

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
    canny_sigma = Parameter(3.0, 'Sigma for the canny algorithm')
    canny_high_threshold = Parameter(0.04, 'High threshold for the canny algorithm')
    canny_low_threshold = Parameter(0.01, 'High threshold for the canny algorithm')

    # Recipe Products
    frame = Product(DataFrameType)
    positions = Product(ArrayType)
    DTU = Product(ArrayType)
    IPA = Product(float)
    param_median_filter_size = Product(float)
    param_canny_high_threshold = Product(float)
    param_canny_low_threshold = Product(float)

    def run(self, rinput):

        logger = logging.getLogger('numina.recipes.emir')

        logger.info('starting processing for bars detection')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            ipa = hdr['IPA']
            xdtu = hdr['XDTU']
            ydtu = hdr['YDTU']
            dtur = get_dtur_from_header(hdr)

        except KeyError as error:
            logger.error(error)
            raise RecipeError(error)

        logger.debug('finding bars')

        arr = hdulist[0].data

        # Median filter
        logger.debug('median filtering')
        mfilter_size = rinput.median_filter_size

        arr_median = median_filter(arr, size=mfilter_size)

        # Image is mapped between 0 and 1
        # for the full range [0: 2**16]
        logger.debug('image scaling to 0-1')
        arr_grey = normalize_raw(arr_median)

        # Find borders
        logger.debug('find borders')
        canny_sigma = rinput.canny_sigma
        # These threshols corespond roughly with
        # value x (2**16 - 1)
        high_threshold = rinput.canny_high_threshold
        low_threshold = rinput.canny_low_threshold

        edges = canny(arr_grey, sigma=canny_sigma,
                      high_threshold=high_threshold,
                      low_threshold=low_threshold)

        # Number or rows used
        # These other parameters cab be tuned also
        total = 5
        maxdist = 1.0
        bstart = 100
        bend = 1900

        positions = []
        nt = total // 2

        xfac = dtur[0] / PIXSCALE
        yfac = -dtur[1] / PIXSCALE

        vec = [yfac, xfac]
        logger.debug('DTU shift is %s', vec)

        # Based on the 'edges image'
        # and the table of approx positions of the slits
        barstab = rinput.bars_nominal_positions

        # Currently, we only use fields 0 and 2
        # of the nominal positions file

        for coords in barstab:
            lbarid = int(coords[0])
            rbarid = lbarid + 55
            ref_y_coor = coords[2] + vec[1]

            logger.debug('looking for bars with ids %i - %i', lbarid, rbarid)
            logger.debug('reference y position is Y %7.2f', ref_y_coor)
            # Find the position of each bar
            bpos = find_position(edges, ref_y_coor, bstart, bend, total, maxdist)

            # If no bar is found, append and empty token
            if bpos is None:
                logger.debug('bar %d not found', lbarid)
                logger.debug('bar %d not found', rbarid)
                thisres1 = (lbarid, 0, 0, 0, 1)
                thisres2 = (rbarid, 0, 0, 0, 1)
            else:
                prow, c1, c2 = bpos
                logger.debug('bars found between %7.2f - %7.2f', c1, c2)
                # Compute FWHM of the collapsed profile

                cslit = arr_grey[prow-nt:prow+nt+1,:]
                pslit = cslit.mean(axis=0)

                # Add 1 to return FITS coordinates
                epos, epos_f, error = locate_bar_l(pslit, c1)
                thisres1 = lbarid, prow + 1, epos + 1, epos_f + 1, error


                epos, epos_f, error = locate_bar_r(pslit, c2)
                thisres2 = rbarid, prow + 1, epos + 1, epos_f + 1, error

            positions.append(thisres1)
            positions.append(thisres2)

        logger.debug('end finding bars')
        result = self.create_result(frame=hdulist,
                                    positions=positions,
                                    DTU=[xdtu, ydtu],
                                    IPA=ipa,
                                    param_median_filter_size=rinput.median_filter_size,
                                    param_canny_high_threshold=rinput.canny_high_threshold,
                                    param_canny_low_threshold=rinput.canny_low_threshold
                                    )
        return result
