#
# Copyright 2008-2018 Universidad Complutense de Madrid
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

from astropy.io import fits
from copy import deepcopy
import logging
import numpy as np

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximplotxy import ximplotxy
from numina.array.stats import summary
from numina.array.wavecalib.fix_pix_borders import find_pix_borders
from numina.array.wavecalib.crosscorrelation import convolve_comb_lines
from numina.array.wavecalib.crosscorrelation import periodic_corr1d

from emirdrp.processing.wavecal.median_slitlets_rectified \
    import median_slitlets_rectified

from emirdrp.core import EMIR_NBARS


def refine_rectwv_coeff(input_image, rectwv_coeff, ohlines, debugplot=0):
    """Refine RectWaveCoeff object using OH sky lines

    Parameters
    ----------
    input_image : HDUList object
        Input 2D image.
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    ohlines : numpy array
        2D numpy array with the contents of the file with the table of
        expected OH sky lines.
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed. The valid codes are defined in
        numina.array.display.pause_debugplot.

    Returns
    -------
    refined_rectwv_coeff_global : RectWaveCoeff instance
        Refined rectification and wavelength calibration coefficients
        for the particular CSU configuration, using the global image
        offset as a correction for CRVAL1.
    refined_rectwv_coeff_individual : RectWaveCoeff instance
        Refined rectification and wavelength calibration coefficients
        for the particular CSU configuration, using the individual
        offset of each slitlet as a correction for CRVAL1.
    expected_oh_image : HDUList object
        Output 2D image with the expected OH sky spectrum.

    """

    logger = logging.getLogger(__name__)

    # initialize output
    refined_rectwv_coeff_global = deepcopy(rectwv_coeff)
    refined_rectwv_coeff_individual = deepcopy(rectwv_coeff)

    # compute median spectrum and normalize it
    sp_median = median_slitlets_rectified(input_image, mode=2)[0].data
    sp_median /= sp_median.max()

    # determine minimum and maximum useful wavelength
    jmin, jmax = find_pix_borders(sp_median, 0)
    main_header = input_image[0].header
    naxis1 = main_header['naxis1']
    naxis2 = main_header['naxis2']
    crpix1 = main_header['crpix1']
    crval1 = main_header['crval1']
    cdelt1 = main_header['cdelt1']
    crvaln = crval1 + (naxis1 - crpix1) * cdelt1
    xwave = crval1 + (np.arange(naxis1) + 1.0 - crpix1) * cdelt1
    wave_min = crval1 + (jmin + 1 - crpix1) * cdelt1
    wave_max = crval1 + (jmax + 1 - crpix1) * cdelt1

    # extract subset of OH lines within current wavelength range
    lok1 = ohlines[:, 1] >= wave_min
    lok2 = ohlines[:, 0] <= wave_max
    ohlines_reference = ohlines[lok1*lok2]

    # define wavelength and flux as separate arrays
    ohlines_reference_wave = np.concatenate(
        (ohlines_reference[:, 1], ohlines_reference[:, 0]))
    ohlines_reference_flux = np.concatenate(
        (ohlines_reference[:, 2], ohlines_reference[:, 2]))
    ohlines_reference_flux /= ohlines_reference_flux.max()

    # convolve location of OH lines to generate expected spectrum
    xwave_reference, sp_oh_reference = convolve_comb_lines(
        ohlines_reference_wave, ohlines_reference_flux, 2.0,
        crpix1, crval1, cdelt1, naxis1
    )
    sp_oh_reference /= sp_oh_reference.max()

    # generate expected_oh_image
    image2d_expected_oh = np.tile(sp_oh_reference, (naxis2, 1))
    hdu = fits.PrimaryHDU(data=image2d_expected_oh, header=main_header)
    expected_oh_image = fits.HDUList([hdu])

    if abs(debugplot) % 10 != 0:
        ax = ximplotxy(xwave, sp_median, 'C0-',
                       xlabel='Wavelength (Angstroms, in vacuum)',
                       ylabel='Normalized number of counts',
                       label='observed spectrum', show=False)
        # overplot reference OH lines
        ax.stem(ohlines_reference_wave, ohlines_reference_flux, 'C4-',
                markerfmt=' ', basefmt='C4-', label='tabulated OH lines')
        # overplot convolved reference OH lines
        ax.plot(xwave_reference, sp_oh_reference, 'C1-',
                label='expected spectrum')

        ax.legend()
        pause_debugplot(debugplot=debugplot, pltshow=True)

    # compute global offset through periodic correlation
    logger.info('Computing global offset')
    offset = periodic_corr1d(
        sp_reference=sp_oh_reference,
        sp_offset=sp_median,
        fminmax=None,
        debugplot=debugplot
    )
    logger.info('Global offset: {} pixels'.format(offset))

    # apply computed offset to obtain refined_rectwv_coeff_global
    missing_slitlets = rectwv_coeff.missing_slitlets
    for islitlet in range(1, EMIR_NBARS + 1):
        if islitlet not in missing_slitlets:
            i = islitlet - 1
            ddum = refined_rectwv_coeff_global.contents[i]
            ddum['wpoly_coeff'][0] -= offset*cdelt1


    # compute individual offset for each slitlet
    median_image = median_slitlets_rectified(input_image, mode=1)
    offset_array = np.zeros(EMIR_NBARS)
    xplot = []
    yplot = []
    for islitlet in range(1, EMIR_NBARS + 1):
        if islitlet not in missing_slitlets:
            i = islitlet - 1
            sp_median = median_image[0].data[i, :]
            sp_median /= sp_median.max()
            offset_array[i] = periodic_corr1d(
                sp_reference=sp_oh_reference,
                sp_offset=median_image[0].data[i, :],
                fminmax=None,
                debugplot=0
            )
            ddum = refined_rectwv_coeff_individual.contents[i]
            ddum['wpoly_coeff'][0] -= offset_array[i]*cdelt1
            xplot.append(islitlet)
            yplot.append(offset_array[i])

    stat_summary = summary(yplot)
    logger.info('Individual offsets (pixels):')
    logger.info('- npoints...:{}'.format(stat_summary['npoints']))
    logger.info('- mean......:{}'.format(stat_summary['mean']))
    logger.info('- median....:{}'.format(stat_summary['median']))
    logger.info('- std.......:{}'.format(stat_summary['std']))
    logger.info('- robust_std:{}'.format(stat_summary['robust_std']))
    if abs(debugplot) % 10 != 0:
        ximplotxy(xplot, yplot,
                  linestyle='', marker='o', color='C0',
                  xlabel='slitlet number',
                  ylabel='offset (pixels)',
                  debugplot=debugplot)

    # return result
    return refined_rectwv_coeff_global, \
           refined_rectwv_coeff_individual,\
           expected_oh_image
