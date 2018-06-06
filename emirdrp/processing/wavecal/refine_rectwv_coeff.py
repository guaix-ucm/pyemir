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


def refine_rectwv_coeff(input_image, rectwv_coeff, catlines, mode, debugplot=0):
    """Refine RectWaveCoeff object using a catalogue of lines

    Parameters
    ----------
    input_image : HDUList object
        Input 2D image.
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    catlines : numpy array
        2D numpy array with the contents of the file with the table of
        expected catalogued lines.
    mode : int
        Integer (from 0 to 2), indicating the type of refinement:
        1 : apply the same global offset to all the slitlets
        2 : apply individual offset to each individual slitlet
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed. The valid codes are defined in
        numina.array.display.pause_debugplot.

    Returns
    -------
    refined_rectwv_coeff : RectWaveCoeff instance
        Refined rectification and wavelength calibration coefficients
        for the particular CSU configuration.
    expected_cat_image : HDUList object
        Output 2D image with the expected catalogue lines.

    """

    logger = logging.getLogger(__name__)

    # protections
    if mode not in [1, 2]:
        raise ValueError('Invalid mode={}'.format(mode))

    # initialize output
    refined_rectwv_coeff = deepcopy(rectwv_coeff)

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
    xwave = crval1 + (np.arange(naxis1) + 1.0 - crpix1) * cdelt1
    wave_min = crval1 + (jmin + 1 - crpix1) * cdelt1
    wave_max = crval1 + (jmax + 1 - crpix1) * cdelt1

    # extract subset of catalogue lines within current wavelength range
    lok1 = catlines[:, 1] >= wave_min
    lok2 = catlines[:, 0] <= wave_max
    catlines_reference = catlines[lok1*lok2]

    # define wavelength and flux as separate arrays
    catlines_reference_wave = np.concatenate(
        (catlines_reference[:, 1], catlines_reference[:, 0]))
    catlines_reference_flux = np.concatenate(
        (catlines_reference[:, 2], catlines_reference[:, 2]))
    catlines_reference_flux /= catlines_reference_flux.max()

    # convolve location of catalogue lines to generate expected spectrum
    xwave_reference, sp_reference = convolve_comb_lines(
        catlines_reference_wave, catlines_reference_flux, 2.0,
        crpix1, crval1, cdelt1, naxis1
    )
    sp_reference /= sp_reference.max()

    # generate expected_oh_image
    image2d_expected_lines = np.tile(sp_reference, (naxis2, 1))
    hdu = fits.PrimaryHDU(data=image2d_expected_lines, header=main_header)
    expected_cat_image = fits.HDUList([hdu])

    if abs(debugplot) % 10 != 0:
        ax = ximplotxy(xwave, sp_median, 'C0-',
                       xlabel='Wavelength (Angstroms, in vacuum)',
                       ylabel='Normalized number of counts',
                       label='observed spectrum', show=False)
        # overplot reference catalogue lines
        ax.stem(catlines_reference_wave, catlines_reference_flux, 'C4-',
                markerfmt=' ', basefmt='C4-', label='tabulated lines')
        # overplot convolved reference lines
        ax.plot(xwave_reference, sp_reference, 'C1-',
                label='expected spectrum')

        ax.legend()
        pause_debugplot(debugplot=debugplot, pltshow=True)

    # compute global offset through periodic correlation
    logger.info('Computing global offset')
    offset = periodic_corr1d(
        sp_reference=sp_reference,
        sp_offset=sp_median,
        fminmax=None,
        debugplot=debugplot
    )
    logger.info('Global offset: {} pixels'.format(offset))

    missing_slitlets = rectwv_coeff.missing_slitlets

    if mode == 1:
        # apply computed offset to obtain refined_rectwv_coeff_global
        for islitlet in range(1, EMIR_NBARS + 1):
            if islitlet not in missing_slitlets:
                i = islitlet - 1
                dumdict = refined_rectwv_coeff.contents[i]
                dumdict['wpoly_coeff'][0] -= offset*cdelt1

    elif mode == 2:
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
                    sp_reference=sp_reference,
                    sp_offset=median_image[0].data[i, :],
                    fminmax=None,
                    debugplot=0
                )
                dumdict = refined_rectwv_coeff.contents[i]
                dumdict['wpoly_coeff'][0] -= offset_array[i]*cdelt1
                xplot.append(islitlet)
                yplot.append(offset_array[i])

        stat_summary = summary(yplot)
        logger.info('Individual offsets (pixels):')
        logger.info('- npoints....: {}'.format(stat_summary['npoints']))
        logger.info('- mean.......: {}'.format(stat_summary['mean']))
        logger.info('- median.....: {}'.format(stat_summary['median']))
        logger.info('- std........: {}'.format(stat_summary['std']))
        logger.info('- robust_std.: {}'.format(stat_summary['robust_std']))
        if abs(debugplot) % 10 != 0:
            ximplotxy(xplot, yplot,
                      linestyle='', marker='o', color='C0',
                      xlabel='slitlet number',
                      ylabel='offset (pixels)',
                      debugplot=debugplot)

    else:

        raise ValueError('Invalide mode={}'.format(mode))

    # return result
    return refined_rectwv_coeff, \
           expected_cat_image
