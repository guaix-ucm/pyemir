#
# Copyright 2015-2018 Universidad Complutense de Madrid
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
from numina.array.utils import coor_to_pix_1d
from numina.core import Requirement, Product, Parameter, RecipeError
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from scipy.ndimage import convolve1d
from scipy.ndimage.filters import median_filter

import emirdrp.datamodel as datamodel
import emirdrp.instrument.distortions as dist
from emirdrp.core import EmirRecipe, EMIR_PIXSCALE, EMIR_NBARS, EMIR_RON
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.products import DataFrameType, NominalPositions
from emirdrp.recipes.aiv.bardetect import char_bar_peak_l, char_bar_peak_r
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement


class BarDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    bars_nominal_positions = Requirement(NominalPositions,
                                         'Nominal positions of the bars'
                                         )
    median_filter_size = Parameter(5, 'Size of the median box')
    average_box_row_size = Parameter(7, 'Number of rows to average for fine centering (odd)')
    average_box_col_size = Parameter(21, 'Number of columns to extract for fine centering (odd)')
    fit_peak_npoints = Parameter(3, 'Number of points to use for fitting the peak (odd)')

    # Recipe Products
    frame = Product(DataFrameType)
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

        self.save_intermediate_img(hdulist, 'reduced_image.fits')

        try:
            rotang = hdr['ROTANG']
            tsutc1 = hdr['TSUTC1']
            dtub, dtur = datamodel.get_dtur_from_header(hdr)
            csupos = datamodel.get_csup_from_header(hdr)
            if len(csupos) != 2 * EMIR_NBARS:
                raise RecipeError('Number of CSUPOS != 2 * NBARS')
            csusens = datamodel.get_cs_from_header(hdr)

        except KeyError as error:
            self.logger.error(error)
            raise RecipeError(error)

        self.logger.debug('start finding bars')
        allpos, slits = self.find_bars(hdulist, rinput, csupos, dtur)
        self.logger.debug('end finding bars')

        if self.intermediate_results:
            with open('ds9.reg', 'w') as ds9reg:
                self.to_ds9_reg(ds9reg, slits)

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

    def median_filtering(self, hdulist, mfilter_size):

        # Processed array
        arr = hdulist[0].data

        self.logger.debug('median filtering 1')
        self.logger.debug('median filtering X, %d columns', mfilter_size)
        arr_median = median_filter(arr, size=(1, mfilter_size))
        self.logger.debug('median filtering X, %d rows', mfilter_size)
        arr_median = median_filter(arr_median, size=(mfilter_size, 1))
        self.save_intermediate_array(arr_median, 'median_image.fits')

        # Median filter of processed array (two times) in the other direction
        # for Y coordinates
        self.logger.debug('median filtering 2')
        self.logger.debug('median filtering Y, %d rows', mfilter_size)
        arr_median_alt = median_filter(arr, size=(mfilter_size, 1))
        self.logger.debug('median filtering Y, %d columns', mfilter_size)
        arr_median_alt = median_filter(arr_median_alt, size=(1, mfilter_size))
        self.save_intermediate_array(arr_median_alt, 'median_image_alt.fits')

        return arr_median, arr_median_alt

    def find_bars(self, hdulist, rinput, csupos, dtur):

        self.logger.debug('filtering image')
        # Processed array
        arr_median, arr_median_alt = self.median_filtering(hdulist,
                                                           rinput.median_filter_size)

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
        self.logger.debug('ignoring columns outside %d - %d', bstart, bend-1)

        # extract a region to average
        wy = (rinput.average_box_row_size // 2)
        wx = (rinput.average_box_col_size // 2)
        self.logger.debug('extraction window is %d rows, %d cols', 2*wy+1, 2*wx+1)
        # Fit the peak with these points
        wfit = 2 * (rinput.fit_peak_npoints // 2) + 1
        self.logger.debug('fit with %d points', wfit)

        # Minimum threshold
        threshold = 5 * EMIR_RON
        # Savitsky and Golay (1964) filter to compute the X derivative
        # scipy >= xx has a savgol_filter function
        # for compatibility we do it manually

        allpos = {}
        slits = numpy.zeros((EMIR_NBARS, 8), dtype='float')

        self.logger.info('find peaks in derivative image')
        for ks in [3, 5, 7, 9]:
            self.logger.debug('kernel size is %d', ks)
            # S and G kernel for derivative
            kw = ks * (ks*ks-1) / 12.0
            coeffs_are = -numpy.arange((1-ks)//2, (ks-1)//2 + 1) / kw
            self.logger.debug('kernel weights are %s', coeffs_are)

            self.logger.debug('derive image in X direction')
            arr_deriv = convolve1d(arr_median, coeffs_are, axis=-1)
            self.save_intermediate_array(arr_deriv, 'deriv_image_k%d.fits' % ks)
            # Axis 0 is
            #
            # self.logger.debug('derive image in Y direction (with kernel=3)')
            # arr_deriv_alt = convolve1d(arr_median_alt, ypos3_kernel, axis=0)

            positions = []
            self.logger.info('using bar parameters')
            for idx in range(EMIR_NBARS):
                params_l = barstab[idx]
                params_r = barstab[idx + EMIR_NBARS]
                lbarid = int(params_l[0])

                # CSUPOS for this bar
                rbarid = lbarid + EMIR_NBARS
                current_csupos_l = csupos[lbarid - 1]
                current_csupos_r = csupos[rbarid - 1]
                self.logger.debug('CSUPOS for bar %d is %f', lbarid, current_csupos_l)
                self.logger.debug('CSUPOS for bar %d is %f', rbarid, current_csupos_r)

                ref_y_coor_virt = params_l[1] # Do I need to add vec[1]?
                ref_x_l_coor_virt = params_l[3] + current_csupos_l * params_l[2]
                ref_x_r_coor_virt = params_r[3] + current_csupos_r * params_r[2]
                # Transform to REAL..
                ref_x_l_coor, ref_y_l_coor = dist.exvp(ref_x_l_coor_virt, ref_y_coor_virt)
                ref_x_r_coor, ref_y_r_coor = dist.exvp(ref_x_r_coor_virt, ref_y_coor_virt)
                # FIXME: check if DTU has to be applied
                # ref_y_coor = ref_y_coor + vec[1]
                prow = coor_to_pix_1d(ref_y_l_coor) - 1
                fits_row = prow + 1 # FITS pixel index

                self.logger.debug('looking for bars with ids %d - %d', lbarid, rbarid)
                self.logger.debug('ref Y virtual position is %7.2f', ref_y_coor_virt)
                self.logger.debug('ref X virtual positions are %7.2f %7.2f', ref_x_l_coor_virt, ref_x_r_coor_virt)
                self.logger.debug('ref X positions are %7.2f %7.2f', ref_x_l_coor, ref_x_r_coor)
                self.logger.debug('ref Y positions are %7.2f %7.2f', ref_y_l_coor, ref_y_r_coor)
                # if ref_y_coor is outlimits, skip this bar
                # ref_y_coor is in FITS format
                if (ref_y_l_coor >= 2047 + 16) or (ref_y_l_coor <= 1 - 16):
                    self.logger.debug('reference y position is outlimits, skipping')
                    positions.append([lbarid, fits_row, fits_row, fits_row, 1, 1, 0, 3])
                    positions.append([rbarid, fits_row, fits_row, fits_row, 1, 1, 0, 3])
                    continue

                # minimal width of the slit
                minwidth = 0.9
                if abs(ref_x_l_coor_virt - ref_x_r_coor_virt) < minwidth:
                    self.logger.debug('slit is less than %d virt pixels, skipping', minwidth)
                    positions.append([lbarid, fits_row, fits_row, fits_row, 1, 1, 0, 3])
                    positions.append([rbarid, fits_row, fits_row, fits_row, 1, 1, 0, 3])
                    continue

                # Left bar
                # Dont add +1 to virtual pixels
                self.logger.debug('measure left border (%d)', lbarid)
                regionw = 10
                bstart1 = coor_to_pix_1d(ref_x_l_coor - regionw)
                bend1 = coor_to_pix_1d(ref_x_l_coor + regionw) + 1
                centery, centery_virt, xpos1, xpos1_virt, fwhm, st = char_bar_peak_l(arr_deriv,
                                                                     prow, bstart1, bend1, threshold,
                                                                     wx=wx, wy=wy, wfit=wfit)

                insert1 = [lbarid, centery + 1, centery_virt, fits_row, xpos1 + 1, xpos1_virt, fwhm, st]
                positions.append(insert1)

                # Right bar
                # Dont add +1 to virtual pixels
                self.logger.debug('measure rigth border (%d)', rbarid)
                bstart2 = coor_to_pix_1d(ref_x_r_coor - regionw)
                bend2 = coor_to_pix_1d(ref_x_r_coor + regionw) + 1
                centery, centery_virt, xpos2, xpos2_virt, fwhm, st = char_bar_peak_r(arr_deriv, prow, bstart2, bend2,
                                                                                     threshold, wx=wx, wy=wy, wfit=wfit)
                # This centery/centery_virt should be equal to ref_y_coor_virt
                insert2 = [rbarid, centery + 1, centery_virt, fits_row, xpos2 + 1, xpos2_virt, fwhm, st]
                positions.append(insert2)

                # FIXME: hardcoded value
                y1_virt = ref_y_coor_virt - 16.242
                y2_virt = ref_y_coor_virt + 16.242
                _, y1 = dist.exvp(xpos1_virt + 1, y1_virt + 1)
                _, y2 = dist.exvp(xpos2_virt + 1, y2_virt + 1)

                # Update positions

                msg = 'bar %d, centroid-y %9.4f centroid-y virt %9.4f, ' \
                      'row %d, x-pos %9.4f x-pos virt %9.4f, FWHM %6.3f, status %d'
                self.logger.debug(msg, *positions[-2])
                self.logger.debug(msg, *positions[-1])

                if ks == 5:
                    slits[lbarid - 1] = numpy.array([xpos1, y2, xpos2, y2, xpos2, y1, xpos1, y1])
                    # FITS coordinates
                    slits[lbarid - 1] += 1.0
                    self.logger.debug('inserting bars %d-%d into "slits"', lbarid, rbarid)

            allpos[ks] = numpy.asarray(positions, dtype='float') # GCS doesn't like lists of lists

        return allpos, slits

    def to_ds9_reg(self, ds9reg, slits):
        """Transform fiber traces to ds9-region format.

        Parameters
        ----------
        ds9reg : BinaryIO
            Handle to output file name in ds9-region format.
        """

        # open output file and insert header

        ds9reg.write('# Region file format: DS9 version 4.1\n')
        ds9reg.write(
            'global color=green dashlist=8 3 width=1 font="helvetica 10 '
            'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 '
            'move=1 delete=1 include=1 source=1\n')
        ds9reg.write('physical\n')

        for idx, slit in enumerate(slits, 1):
            xpos1, y2, xpos2, y2, xpos2, y1, xpos1, y1 = slit
            xc = 0.5 * (xpos1 + xpos2) + 1
            yc = 0.5 * (y1 + y2) + 1
            xd = (xpos2 - xpos1)
            yd = (y2 - y1)
            ds9reg.write('box({0},{1},{2},{3},0)\n'.format(xc, yc, xd, yd))
            ds9reg.write('# text({0},{1}) color=red text={{{2}}}\n'.format(xpos1 - 5, yc, idx))
            ds9reg.write('# text({0},{1}) color=red text={{{2}}}\n'.format(xpos2 + 5, yc, idx + EMIR_NBARS))
