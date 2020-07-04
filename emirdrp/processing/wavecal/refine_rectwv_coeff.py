#
# Copyright 2008-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import division
from __future__ import print_function

from copy import deepcopy
import logging
import numpy as np

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximplotxy import ximplotxy
from numina.array.display.matplotlib_qt import plt
from numina.array.display.matplotlib_qt import set_window_geometry
from numina.array.stats import summary
from numina.array.wavecalib.check_wlcalib import check_wlcalib_sp
from numina.array.wavecalib.crosscorrelation import convolve_comb_lines
from numina.array.wavecalib.crosscorrelation import periodic_corr1d
from numina.array.wavecalib.fix_pix_borders import find_pix_borders
from numina.frame.utils import copy_img

import emirdrp.datamodel as datamodel
from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.processing.wavecal.median_slitlets_rectified \
    import median_slitlets_rectified
from emirdrp.processing.wavecal.set_wv_parameters import set_wv_parameters

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NBARS


def refine_rectwv_coeff(input_image, rectwv_coeff,
                        catlines_all_wave,
                        catlines_all_flux,
                        refine_wavecalib_mode,
                        list_useful_slitlets,
                        save_intermediate_results=False,
                        debugplot=0):
    """Refine RectWaveCoeff object using a catalogue of lines

    One and only one among refine_with_oh_lines_mode and
    refine_with_arc_lines must be different from zero.

    Parameters
    ----------
    input_image : HDUList object
        Input 2D image.
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    catlines_all_wave : numpy array
        Array with wavelengths
    catlines_all_flux : numpy array
        Array with fluxes.
    refine_wavecalib_mode : int
        Integer, indicating the type of refinement:
        0 : no refinement
        1 : apply the same global offset to all the slitlets (using ARC lines)
        2 : apply individual offset to each slitlet (using ARC lines)
        11 : apply the same global offset to all the slitlets (using OH lines)
        12 : apply individual offset to each slitlet (using OH lines)
    list_useful_slitlets : list of integers
        List of useful slitlets.
    save_intermediate_results : bool
        If True, save plots in PDF files
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
        Output 2D image (rectified and wavelength calibrated) with
        the expected catalogue lines.

    """

    logger = logging.getLogger(__name__)

    if save_intermediate_results:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages('crosscorrelation.pdf')
    else:
        pdf = None

    # image header
    main_header = input_image[0].header
    mecs_header = datamodel.get_mecs_header(input_image)
    filter_name = main_header['filter']
    grism_name = main_header['grism']

    # initialize output
    refined_rectwv_coeff = deepcopy(rectwv_coeff)

    logger.info('Computing median spectrum')
    # compute median spectrum and normalize it
    sp_median = median_slitlets_rectified(
        input_image,
        mode=2,
        list_useful_slitlets=list_useful_slitlets
    )[0].data
    sp_median /= sp_median.max()

    # determine minimum and maximum useful wavelength
    jmin, jmax = find_pix_borders(sp_median, 0)
    naxis1 = main_header['naxis1']
    naxis2 = main_header['naxis2']
    crpix1 = main_header['crpix1']
    crval1 = main_header['crval1']
    cdelt1 = main_header['cdelt1']
    xwave = crval1 + (np.arange(naxis1) + 1.0 - crpix1) * cdelt1
    if grism_name == 'LR':
        wv_parameters = set_wv_parameters(filter_name, grism_name)
        wave_min = wv_parameters['wvmin_useful']
        wave_max = wv_parameters['wvmax_useful']
    else:
        wave_min = crval1 + (jmin + 1 - crpix1) * cdelt1
        wave_max = crval1 + (jmax + 1 - crpix1) * cdelt1
    logger.info('Setting wave_min to {}'.format(wave_min))
    logger.info('Setting wave_max to {}'.format(wave_max))

    # extract subset of catalogue lines within current wavelength range
    lok1 = catlines_all_wave >= wave_min
    lok2 = catlines_all_wave <= wave_max
    catlines_reference_wave = catlines_all_wave[lok1*lok2]
    catlines_reference_flux = catlines_all_flux[lok1*lok2]
    catlines_reference_flux /= catlines_reference_flux.max()

    # estimate sigma to broaden catalogue lines
    csu_config = CsuConfiguration.define_from_header(mecs_header)
    # segregate slitlets
    list_not_useful_slitlets = [i for i in list(range(1, EMIR_NBARS + 1))
                                if i not in list_useful_slitlets]
    logger.info('list of useful slitlets: {}'.format(
        list_useful_slitlets))
    logger.info('list of unusable slitlets: {}'.format(
        list_not_useful_slitlets))
    tempwidths = np.array([csu_config.csu_bar_slit_width(islitlet)
                           for islitlet in list_useful_slitlets])
    widths_summary = summary(tempwidths)
    logger.info('Statistics of useful slitlet widths (mm):')
    logger.info('- npoints....: {0:d}'.format(widths_summary['npoints']))
    logger.info('- mean.......: {0:7.3f}'.format(widths_summary['mean']))
    logger.info('- median.....: {0:7.3f}'.format(widths_summary['median']))
    logger.info('- std........: {0:7.3f}'.format(widths_summary['std']))
    logger.info('- robust_std.: {0:7.3f}'.format(widths_summary['robust_std']))
    # empirical transformation of slit width (mm) to pixels
    sigma_broadening = cdelt1 * widths_summary['median']

    # convolve location of catalogue lines to generate expected spectrum
    xwave_reference, sp_reference = convolve_comb_lines(
        catlines_reference_wave, catlines_reference_flux, sigma_broadening,
        crpix1, crval1, cdelt1, naxis1
    )
    sp_reference /= sp_reference.max()

    # generate image2d with expected lines
    if save_intermediate_results:
        image2d_expected_lines = np.tile(sp_reference, (naxis2, 1))
        expected_cat_image = copy_img(input_image)
        expected_cat_image[0].data = image2d_expected_lines
    else:
        expected_cat_image = None

    if (abs(debugplot) % 10 != 0) or (pdf is not None):
        ax = ximplotxy(xwave, sp_median, 'C1-',
                       xlabel='Wavelength (Angstroms, in vacuum)',
                       ylabel='Normalized number of counts',
                       title='Median spectrum',
                       label='observed spectrum', show=False)
        # overplot reference catalogue lines
        ax.stem(catlines_reference_wave, catlines_reference_flux, 'C4-',
                markerfmt=' ', basefmt='C4-', label='tabulated lines')
        # overplot convolved reference lines
        ax.plot(xwave_reference, sp_reference, 'C0-',
                label='expected spectrum')
        ax.legend()
        if pdf is not None:
            pdf.savefig()
        else:
            pause_debugplot(debugplot=debugplot, pltshow=True)

    # compute baseline signal in sp_median
    baseline = np.percentile(sp_median[sp_median > 0], q=10)
    if (abs(debugplot) % 10 != 0) or (pdf is not None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.hist(sp_median, bins=1000, log=True)
        ax.set_xlabel('Normalized number of counts')
        ax.set_ylabel('Number of pixels')
        ax.set_title('Median spectrum')
        ax.axvline(float(baseline), linestyle='--', color='grey')
        if pdf is not None:
            pdf.savefig()
        else:
            geometry = (0, 0, 640, 480)
            set_window_geometry(geometry)
            plt.show()
    # subtract baseline to sp_median (only pixels with signal above zero)
    lok = np.where(sp_median > 0)
    sp_median[lok] -= baseline

    # compute global offset through periodic correlation
    logger.info('Computing global offset')
    global_offset, fpeak = periodic_corr1d(
        sp_reference=sp_reference,
        sp_offset=sp_median,
        fminmax=None,
        naround_zero=50,
        plottitle='Median spectrum (cross-correlation)',
        pdf=pdf,
        debugplot=debugplot
    )
    logger.info('Global offset: {} pixels'.format(-global_offset))

    missing_slitlets = rectwv_coeff.missing_slitlets

    if refine_wavecalib_mode in [1, 11]:
        # apply computed offset to obtain refined_rectwv_coeff_global
        for islitlet in range(1, EMIR_NBARS + 1):
            if islitlet not in missing_slitlets:
                i = islitlet - 1
                dumdict = refined_rectwv_coeff.contents[i]
                dumdict['wpoly_coeff'][0] -= global_offset*cdelt1

    elif refine_wavecalib_mode in [2, 12]:
        # compute individual offset for each slitlet
        logger.info('Computing individual offsets')
        median_55sp = median_slitlets_rectified(input_image, mode=1)
        offset_array = np.zeros(EMIR_NBARS)
        xplot = []
        yplot = []
        xplot_skipped = []
        yplot_skipped = []
        cout = '0'
        for islitlet in range(1, EMIR_NBARS + 1):
            if islitlet in list_useful_slitlets:
                i = islitlet - 1
                sp_median = median_55sp[0].data[i, :]
                lok = np.where(sp_median > 0)
                if np.any(lok):
                    baseline = np.percentile(sp_median[lok], q=10)
                    sp_median[lok] -= baseline
                    sp_median /= sp_median.max()
                    offset_array[i], fpeak = periodic_corr1d(
                        sp_reference=sp_reference,
                        sp_offset=median_55sp[0].data[i, :],
                        fminmax=None,
                        naround_zero=50,
                        plottitle='slitlet #{0} (cross-correlation)'.format(
                            islitlet),
                        pdf=pdf,
                        debugplot=debugplot
                    )
                else:
                    offset_array[i] = 0.0
                dumdict = refined_rectwv_coeff.contents[i]
                dumdict['wpoly_coeff'][0] -= offset_array[i]*cdelt1
                xplot.append(islitlet)
                yplot.append(-offset_array[i])
                # second correction
                wpoly_coeff_refined = check_wlcalib_sp(
                    sp=median_55sp[0].data[i, :],
                    crpix1=crpix1,
                    crval1=crval1-offset_array[i]*cdelt1,
                    cdelt1=cdelt1,
                    wv_master=catlines_reference_wave,
                    coeff_ini=dumdict['wpoly_coeff'],
                    naxis1_ini=EMIR_NAXIS1,
                    title='slitlet #{0} (after applying offset)'.format(
                        islitlet),
                    ylogscale=False,
                    pdf=pdf,
                    debugplot=debugplot
                )
                dumdict['wpoly_coeff'] = wpoly_coeff_refined
                cout += '.'

            else:
                xplot_skipped.append(islitlet)
                yplot_skipped.append(0)
                cout += 'i'

            if islitlet % 10 == 0:
                if cout != 'i':
                    cout = str(islitlet // 10)

            logger.info(cout)

        # show offsets with opposite sign
        stat_summary = summary(np.array(yplot))
        logger.info('Statistics of individual slitlet offsets (pixels):')
        logger.info('- npoints....: {0:d}'.format(stat_summary['npoints']))
        logger.info('- mean.......: {0:7.3f}'.format(stat_summary['mean']))
        logger.info('- median.....: {0:7.3f}'.format(stat_summary['median']))
        logger.info('- std........: {0:7.3f}'.format(stat_summary['std']))
        logger.info('- robust_std.: {0:7.3f}'.format(stat_summary[
                                                        'robust_std']))
        if (abs(debugplot) % 10 != 0) or (pdf is not None):
            ax = ximplotxy(xplot, yplot,
                           linestyle='', marker='o', color='C0',
                           xlabel='slitlet number',
                           ylabel='-offset (pixels) = offset to be applied',
                           title='cross-correlation result',
                           show=False, **{'label': 'individual slitlets'})
            if len(xplot_skipped) > 0:
                ax.plot(xplot_skipped, yplot_skipped, 'mx')
            ax.axhline(-global_offset, linestyle='--', color='C1',
                       label='global offset')
            ax.legend()
            if pdf is not None:
                pdf.savefig()
            else:
                pause_debugplot(debugplot=debugplot, pltshow=True)
    else:
        raise ValueError('Unexpected mode={}'.format(refine_wavecalib_mode))

    # close output PDF file
    if pdf is not None:
        pdf.close()

    # return result
    return refined_rectwv_coeff, expected_cat_image
