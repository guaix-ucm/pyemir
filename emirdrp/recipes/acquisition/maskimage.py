#
# Copyright 2015-2016 Universidad Complutense de Madrid
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

import numpy
from numina.array.utils import wc_to_pix_1d
from numina.core import Requirement, Product, Parameter, RecipeError
from numina.core.products import ArrayType
from numpy.polynomial.polynomial import polyval
from scipy.ndimage import convolve1d
from scipy.ndimage.filters import median_filter

import emirdrp.datamodel as datamodel
import emirdrp.products as prods
import emirdrp.requirements as reqs
from emirdrp.core import EmirRecipe, EMIR_PIXSCALE, EMIR_NBARS, EMIR_RON
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.recipes.aiv.bardetect import char_bar_peak_l, char_bar_peak_r, char_bar_height


class MaskImagingRecipe(EmirRecipe):

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
    average_box_row_size = Parameter(7, 'Number of rows to average for fine centering (odd)')
    average_box_col_size = Parameter(21, 'Number of columns to extract for fine centering (odd)')
    fit_peak_npoints = Parameter(3, 'Number of points to use for fitting the peak (odd)')

    # Recipe Products
    frame = Product(prods.DataFrameType)
    # derivative = Product(DataFrameType)
    slits = Product(ArrayType)
    positions3 = Product(ArrayType)
    positions5 = Product(ArrayType)
    positions7 = Product(ArrayType)
    positions9 = Product(ArrayType)
    DTU = Product(ArrayType)
    ROTANG = Product(float)
    TSUTC1 = Product(float)
    csupos = Product(ArrayType)
    csusens = Product(ArrayType)

    def run(self, rinput):
        self.logger.info('starting processing for bars detection')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            rotang = hdr['ROTANG']
            tsutc1 = hdr['TSUTC1']
            dtub, dtur = datamodel.get_dtur_from_header(hdr)
            csupos = datamodel.get_csup_from_header(hdr)
            csusens = datamodel.get_cs_from_header(hdr)

        except KeyError as error:
            self.logger.error(error)
            raise RecipeError(error)

        self.logger.debug('finding bars')
        # Processed array
        arr = hdulist[0].data

        # Median filter of processed array (two times)
        mfilter_size = rinput.median_filter_size

        self.logger.debug('median filtering X, %d columns', mfilter_size)
        arr_median = median_filter(arr, size=(1, mfilter_size))
        self.logger.debug('median filtering X, %d rows', mfilter_size)
        arr_median = median_filter(arr_median, size=(mfilter_size, 1))

        # Median filter of processed array (two times) in the other direction
        # for Y coordinates
        self.logger.debug('median filtering Y, %d rows', mfilter_size)
        arr_median_alt = median_filter(arr, size=(mfilter_size, 1))
        self.logger.debug('median filtering Y, %d columns', mfilter_size)
        arr_median_alt = median_filter(arr_median_alt, size=(1, mfilter_size))

        xfac = dtur[0] / EMIR_PIXSCALE
        yfac = -dtur[1] / EMIR_PIXSCALE

        vec = [yfac, xfac]
        self.logger.debug('DTU shift is %s', vec)

        # and the table of approx positions of the slits
        barstab = rinput.bars_nominal_positions
        # Currently, we only use fields 0 and 2
        # of the nominal positions file

        # Number or rows used
        # These other parameters cab be tuned also
        bstart = 1
        bend = 2047
        self.logger.debug('ignoring columns outside %d - %d',bstart, bend-1)

        # extract a region to average
        wy = (rinput.average_box_row_size // 2)
        wx = (rinput.average_box_col_size // 2)
        self.logger.debug('extraction window is %d rows, %d cols',2*wy+1, 2*wx+1)
        # Fit the peak with these points
        wfit = 2 * (rinput.fit_peak_npoints // 2) + 1
        self.logger.debug('fit with %d points', wfit)

        # Minimum threshold
        threshold = 5 * EMIR_RON
        # Savitsky and Golay (1964) filter to compute the X derivative
        # scipy >= xx has a savgol_filter function
        # for compatibility we do it manually

        allpos = {}
        ypos3_kernel = None
        slits = numpy.zeros((EMIR_NBARS, 8), dtype='float')

        self.logger.info('start finding bars')
        for ks in [3, 5, 7, 9]:
            self.logger.debug('kernel size is %d', ks)
            # S and G kernel for derivative
            kw = ks * (ks*ks-1) / 12.0
            coeffs_are = -numpy.arange((1-ks)//2, (ks-1)//2 + 1) / kw
            if ks == 3:
                ypos3_kernel = coeffs_are
            self.logger.debug('kernel weights are %s', coeffs_are)

            self.logger.debug('derive image in X direction')
            arr_deriv = convolve1d(arr_median, coeffs_are, axis=-1)
            # Axis 0 is
            #
            self.logger.debug('derive image in Y direction (with kernel=3)')
            arr_deriv_alt = convolve1d(arr_median_alt, ypos3_kernel, axis=0)

            positions = []
            for coords in barstab:
                lbarid = int(coords[0])
                rbarid = lbarid + EMIR_NBARS
                ref_y_coor = coords[1] + vec[1]
                poly_coeffs = coords[2:]
                prow = wc_to_pix_1d(ref_y_coor) - 1
                fits_row = prow + 1 # FITS pixel index

                # A function that returns the center of the bar
                # given its X position
                def center_of_bar(x):
                    # Pixel values are 0-based
                    return polyval(x+1-vec[0], poly_coeffs) + vec[1] - 1

                self.logger.debug('looking for bars with ids %d - %d', lbarid, rbarid)
                self.logger.debug('reference y position is Y %7.2f', ref_y_coor)

                # if ref_y_coor is outlimits, skip this bar
                # ref_y_coor is in FITS format
                if (ref_y_coor >= 2047) or (ref_y_coor <= 1):
                    self.logger.debug('reference y position is outlimits, skipping')
                    positions.append([lbarid, fits_row, fits_row, 1, 0, 3])
                    positions.append([rbarid, fits_row, fits_row, 1, 0, 3])
                    continue

                # Left bar
                self.logger.debug('measure left border (%d)', lbarid)

                centery, xpos, fwhm, st = char_bar_peak_l(arr_deriv, prow, bstart, bend, threshold,
                                                          center_of_bar, wx=wx, wy=wy, wfit=wfit)
                xpos1 = xpos
                positions.append([lbarid, centery+1, fits_row, xpos+1, fwhm, st])

                # Right bar
                self.logger.debug('measure rigth border (%d)', rbarid)
                centery, xpos, fwhm, st = char_bar_peak_r(arr_deriv, prow, bstart, bend, threshold,
                                                          center_of_bar, wx=wx, wy=wy, wfit=wfit)
                positions.append([rbarid, centery+1, fits_row, xpos+1, fwhm, st])
                xpos2 = xpos
                #
                if st == 0:
                    self.logger.debug('measure top-bottom borders')
                    try:
                        y1, y2, statusy = char_bar_height(arr_deriv_alt, xpos1, xpos2, centery, threshold,
                                                          wh=35, wfit=wfit)
                    except Exception as error:
                        self.logger.warning('Error computing height: %s', error)
                        statusy = 44

                    if statusy in [0, 40]:
                        # Main border is detected
                        positions[-1][1] = y2 + 1
                        positions[-2][1] = y2 + 1
                    else:
                        # Update status
                        positions[-1][-1] = 4
                        positions[-2][-1] = 4
                else:
                    self.logger.debug('slit is not complete')
                    y1, y2 = 0, 0

                # Update positions

                self.logger.debug('bar %d centroid-y %9.4f, row %d x-pos %9.4f, FWHM %6.3f, status %d', *positions[-2])
                self.logger.debug('bar %d centroid-y %9.4f, row %d x-pos %9.4f, FWHM %6.3f, status %d', *positions[-1])

                if ks == 5:
                    slits[lbarid - 1] = [xpos1, y2, xpos2, y2, xpos2, y1, xpos1, y1]
                    # FITS coordinates
                    slits[lbarid - 1] += 1.0
                    self.logger.debug('inserting bars %d-%d into "slits"', lbarid, rbarid)

            allpos[ks] = numpy.asarray(positions, dtype='float') # GCS doesn't like lists of lists

        self.logger.debug('end finding bars')
        result = self.create_result(frame=hdulist,
                                    slits=slits,
                                    positions9=allpos[9],
                                    positions7=allpos[7],
                                    positions5=allpos[5],
                                    positions3=allpos[3],
                                    DTU=dtub,
                                    ROTANG=rotang,
                                    TSUTC1=tsutc1,
                                    csupos=csupos,
                                    csusens=csusens,
                                    )
        return result
