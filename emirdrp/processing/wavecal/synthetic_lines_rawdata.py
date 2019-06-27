#
# Copyright 2019 Universidad Complutense de Madrid
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

from __future__ import division
from __future__ import print_function

import logging
import numpy as np
from scipy.stats import norm

from numina.array.distortion import fmap

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS


def synthetic_lines_rawdata(catlines_all_wave,
                            catlines_all_flux,
                            list_useful_slitlets,
                            rectwv_coeff,
                            ycorrection_pix=2.0,
                            sigma_gauss_pix = 1.5):
    """Generate synthetic MOS data with arc/OH lines.

    Parameters
    ----------
    catlines_all_wave : numpy array
        Array with wavelengths of the lines to be overplotted.
    catlines_all_flux : numpy array
        Array with fluxes (arbitrary units) of the lines to be
        overplotted.
    list_useful_slitlets : list of integers
        List with numbers of valid slitlets.
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    ycorrection_pix : float
        Approximate number of pixels from frontier to boundary.
    sigma_gauss_pix : float
        Sigma of the Gaussians to be employed to reproduce the
        emission lines.

    Returns
    -------
    simulated_image : numpy array
        Simulated image with arc/OH lines.

    """

    logger = logging.getLogger(__name__)
    logger.info('generating synthetic MOS image')

    simulated_image = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

    global_integer_offset_x_pix = rectwv_coeff.global_integer_offset_x_pix
    global_integer_offset_y_pix = rectwv_coeff.global_integer_offset_y_pix

    nside_pix = int(sigma_gauss_pix * 5.0 + 0.5)

    cout = '0'
    for islitlet in range(1, EMIR_NBARS + 1):
        if islitlet in list_useful_slitlets:
            cout += '.'
            dumdict = rectwv_coeff.contents[islitlet - 1]
            # check pseudo-longslit with previous slitlet
            center = dumdict['csu_bar_slit_center']
            width = dumdict['csu_bar_slit_width']
            connected_prev = False
            if islitlet > 1:
                dumdict_prev = rectwv_coeff.contents[islitlet - 2]
                center_prev = dumdict_prev['csu_bar_slit_center']
                width_prev = dumdict_prev['csu_bar_slit_width']
                if abs(center - center_prev) < width/2:
                    if abs(width - width_prev) < width/2:
                        connected_prev = True
            # check pseudo-longslit with next slitlet
            connected_next = False
            if islitlet < EMIR_NBARS:
                dumdict_next = rectwv_coeff.contents[islitlet]
                center_next = dumdict_next['csu_bar_slit_center']
                width_next = dumdict_next['csu_bar_slit_width']
                if abs(center - center_next) < width/2:
                    if abs(width - width_next) < width/2:
                        connected_next = True
            # wavelength calibration parameters
            crval1_linear = dumdict['crval1_linear']
            cdelt1_linear = dumdict['cdelt1_linear']
            crvaln_linear = crval1_linear + (EMIR_NAXIS1 - 1) * cdelt1_linear
            # rectification parameters
            bb_ns1_orig = dumdict['bb_ns1_orig']
            ttd_order = dumdict['ttd_order']
            aij = dumdict['ttd_aij']
            bij = dumdict['ttd_bij']
            min_row_rectified = float(dumdict['min_row_rectified'])
            max_row_rectified = float(dumdict['max_row_rectified'])
            wpoly_coeff = dumdict['wpoly_coeff']
            # store coordinates for each line
            x1 = []
            y1 = []
            x2 = []
            y2 = []
            flist = []
            for wave, flux in zip(catlines_all_wave, catlines_all_flux):
                if crval1_linear <= wave <= crvaln_linear:
                    tmp_coeff = np.copy(wpoly_coeff)
                    tmp_coeff[0] -= wave
                    tmp_xroots = np.polynomial.Polynomial(tmp_coeff).roots()
                    for dum in tmp_xroots:
                        if np.isreal(dum):
                            dum = dum.real
                            if 1 <= dum <= EMIR_NAXIS1:
                                x1.append(dum)
                                y1.append(min_row_rectified)
                                x2.append(dum)
                                y2.append(max_row_rectified)
                                flist.append(flux)
            xx1, yy1 = fmap(ttd_order, aij, bij, np.array(x1), np.array(y1))
            xx1 -= global_integer_offset_x_pix
            yy1 += bb_ns1_orig - global_integer_offset_y_pix
            if not connected_prev:
                yy1 += ycorrection_pix
            xx2, yy2 = fmap(ttd_order, aij, bij, np.array(x2), np.array(y2))
            xx2 -= global_integer_offset_x_pix
            yy2 += bb_ns1_orig - global_integer_offset_y_pix
            if not connected_next:
                yy2 -= ycorrection_pix

            for xx1_, xx2_, yy1_, yy2_, flux in zip(xx1, xx2, yy1, yy2, flist):
                slope = (xx2_ - xx1_) / (yy2_ - yy1_)
                iyy1_ = int(yy1_ + 0.5)
                if yy1_ <= float(iyy1_):
                    fracpix1 = (float(iyy1_) - yy1_)
                    iyy1_ -= 1
                else:
                    fracpix1 = 1.0 - (yy1_ - float(iyy1_))
                iyy2_ = int(yy2_ + 0.5)
                if yy2_ <= float(iyy2_):
                    fracpix2 = 1.0 - (float(iyy2_) - yy2_)
                else:
                    fracpix2 = (yy2_ - float(iyy2_))
                    iyy2_ += 1
                for iyy_ in range(iyy1_, iyy2_ + 1):
                    if iyy_ == iyy1_:
                        fracpix = fracpix1
                    elif iyy_ == iyy2_:
                        fracpix = fracpix2
                    else:
                        fracpix = 1.0
                    icenter = iyy_ - 1
                    if 0 <= icenter <= EMIR_NAXIS2 - 1:
                        xx0_ = xx1_ + slope * (iyy_ - yy1_)
                        xcenter = int(xx0_ + 0.5)
                        npix = 2 * nside_pix + 1
                        xpix = xcenter + np.arange(-nside_pix, nside_pix + 1)
                        left_border = xpix - xx0_ - 0.5
                        right_border = left_border + 1.0
                        borders = np.concatenate((left_border,
                                                  right_border))
                        area = norm.cdf(borders, scale=sigma_gauss_pix)
                        for i in range(npix):
                            xdum = xpix[i]
                            if 1 <= xdum <= EMIR_NAXIS1:
                                simulated_image[icenter, xdum - 1] += \
                                    fracpix * flux * (area[i + npix] - area[i])
        else:
            cout += 'i'

        if islitlet % 10 == 0:
            if cout != 'i':
                cout = str(islitlet // 10)

        logger.info(cout)

    return simulated_image
