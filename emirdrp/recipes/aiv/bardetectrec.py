#
# Copyright 2015-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Bar detection recipe for EMIR"""

from __future__ import division

import logging

from numina.array.utils import coor_to_pix_1d
from numina.core import Requirement, Result, Parameter
import numina.exceptions
import numina.types.array as tarray
from numina.processing.combine import basic_processing_with_combination
from scipy.ndimage.filters import median_filter
from skimage.feature import canny

import emirdrp.datamodel as datamodel
from emirdrp.core import EMIR_PIXSCALE
from emirdrp.core.recipe import EmirRecipe
import emirdrp.requirements as reqs
import emirdrp.products as prods
from emirdrp.processing.bardetect import find_position
from emirdrp.processing.bardetect import locate_bar_l, locate_bar_r
from .common import normalize_raw


class BarDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement()

    bars_nominal_positions = Requirement(prods.CoordinateList2DType,
                                         'Nominal positions of the bars'
                                         )
    median_filter_size = Parameter(5, 'Size of the median box')
    canny_sigma = Parameter(3.0, 'Sigma for the canny algorithm')
    canny_high_threshold = Parameter(0.04, 'High threshold for the canny algorithm')
    canny_low_threshold = Parameter(0.01, 'High threshold for the canny algorithm')

    # Recipe Results
    frame = Result(prods.ProcessedImage)
    positions = Result(tarray.ArrayType)
    DTU = Result(tarray.ArrayType)
    ROTANG = Result(float)
    csupos = Result(tarray.ArrayType)
    csusens = Result(tarray.ArrayType)
    param_median_filter_size = Result(float)
    param_canny_high_threshold = Result(float)
    param_canny_low_threshold = Result(float)

    def run(self, rinput):

        logger = logging.getLogger('numina.recipes.emir')

        logger.info('starting processing for bars detection')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            rotang = hdr['ROTANG']
            mecs_hdr = datamodel.get_mecs_header(hdulist)
            dtub, dtur = datamodel.get_dtur_from_header(mecs_hdr)
            csupos = datamodel.get_csup_from_header(mecs_hdr)
            csusens = datamodel.get_cs_from_header(mecs_hdr)
        except KeyError as error:
            logger.error(error)
            raise numina.exceptions.RecipeError(error)

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

        xfac = dtur[0] / EMIR_PIXSCALE
        yfac = -dtur[1] / EMIR_PIXSCALE

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
            prow = coor_to_pix_1d(ref_y_coor) - 1
            fits_row = prow + 1 # FITS pixel index

            logger.debug('looking for bars with ids %d - %d', lbarid, rbarid)
            logger.debug('reference y position is Y %7.2f', ref_y_coor)
            # Find the position of each bar

            bpos = find_position(edges, prow, bstart, bend, total)

            nbars_found = len(bpos)

            # If no bar is found, append and empty token
            if nbars_found == 0:
                logger.debug('bars %d, %d not found at row %d', lbarid, rbarid, fits_row)
                thisres1 = (lbarid, fits_row, 0, 0, 1)
                thisres2 = (rbarid, fits_row, 0, 0, 1)

            elif nbars_found == 2:

                # Order values by increasing X
                centl, centr = sorted(bpos, key=lambda cen: cen[0])
                c1 = centl[0]
                c2 = centr[0]

                logger.debug('bars found  at row %d between %7.2f - %7.2f', fits_row, c1, c2)
                # Compute FWHM of the collapsed profile

                cslit = arr_grey[prow-nt:prow+nt+1,:]
                pslit = cslit.mean(axis=0)

                # Add 1 to return FITS coordinates
                epos, epos_f, error = locate_bar_l(pslit, c1)
                thisres1 = lbarid, fits_row, epos + 1, epos_f + 1, error


                epos, epos_f, error = locate_bar_r(pslit, c2)
                thisres2 = rbarid, fits_row, epos + 1, epos_f + 1, error

            elif nbars_found == 1:
                logger.warning('only 1 edge found  at row %d, not yet implemented', fits_row)
                thisres1 = (lbarid, fits_row, 0, 0, 1)
                thisres2 = (rbarid, fits_row, 0, 0, 1)

            else:
                logger.warning('3 or more edges found  at row %d, not yet implemented', fits_row)
                thisres1 = (lbarid, fits_row, 0, 0, 1)
                thisres2 = (rbarid, fits_row, 0, 0, 1)

            positions.append(thisres1)
            positions.append(thisres2)

        logger.debug('end finding bars')
        result = self.create_result(frame=hdulist,
                                    positions=positions,
                                    DTU=dtub,
                                    ROTANG=rotang,
                                    csupos=csupos,
                                    csusens=csusens,
                                    param_median_filter_size=rinput.median_filter_size,
                                    param_canny_high_threshold=rinput.canny_high_threshold,
                                    param_canny_low_threshold=rinput.canny_low_threshold
                                    )
        return result
