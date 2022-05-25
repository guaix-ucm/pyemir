#
# Copyright 2018-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Rectification and wavelength calibration polynomials from arc image
"""

from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
from datetime import datetime
import logging
import numpy as np
from scipy.signal import medfilt
import sys
from uuid import uuid4

from numina.array.display.logging_from_debugplot import logging_from_debugplot
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.wavecalib.__main__ import read_wv_master_from_array
from numina.array.wavecalib.__main__ import wvcal_spectrum
from numina.array.wavecalib.arccalibration import refine_arccalibration
from numina.frame.utils import copy_img
from numina.tools.arg_file_is_new import arg_file_is_new
import numina.types.qc

import emirdrp.datamodel as datamodel
from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.instrument.dtuconf import DtuConf
from emirdrp.products import RefinedBoundaryModelParam
from emirdrp.products import RectWaveCoeff
from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.select_unrectified_slitlets import \
    select_unrectified_slitlet

from .set_wv_parameters import set_wv_parameters
from .slitlet2darc import Slitlet2dArc

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NBARS


def rectwv_coeff_from_arc_image(reduced_image,
                                bound_param,
                                lines_catalog,
                                args_nbrightlines=None,
                                args_ymargin_bb=2,
                                args_remove_sp_background=True,
                                args_times_sigma_threshold=10,
                                args_order_fmap=2,
                                args_sigma_gaussian_filtering=2,
                                args_margin_npix=50,
                                args_poldeg_initial=3,
                                args_poldeg_refined=5,
                                args_interactive=False,
                                args_threshold_wv=0,
                                args_ylogscale=False,
                                args_pdf=None,
                                args_geometry=(0,0,640,480),
                                debugplot=0):
    """Evaluate rect.+wavecal. coefficients from arc image

    Parameters
    ----------
    reduced_image : HDUList object
        Image with preliminary basic reduction: bpm, bias, dark and
        flatfield.
    bound_param : RefinedBoundaryModelParam instance
        Refined boundary model.
    lines_catalog : Numpy array
        2D numpy array with the contents of the master file with the
        expected arc line wavelengths.
    args_nbrightlines : int
        TBD
    args_ymargin_bb : int
        TBD
    args_remove_sp_background : bool
        TBD
    args_times_sigma_threshold : float
        TBD
    args_order_fmap : int
        TBD
    args_sigma_gaussian_filtering : float
        TBD
    args_margin_npix : int
        TBD
    args_poldeg_initial : int
        TBD
    args_poldeg_refined : int
        TBD
    args_interactive : bool
        TBD
    args_threshold_wv : float
        TBD
    args_ylogscale : bool
        TBD
    args_pdf : TBD
    args_geometry : TBD
    debugplot : int
            Debugging level for messages and plots. For details see
            'numina.array.display.pause_debugplot.py'.

    Returns
    -------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration of the input arc image.
    reduced_55sp : HDUList object
        Image with 55 spectra corresponding to the median spectrum for
        each slitlet, employed to derived the wavelength calibration
        polynomial.

    """

    logger = logging.getLogger(__name__)

    # protections
    if args_interactive and args_pdf is not None:
        logger.error('--interactive and --pdf are incompatible options')
        raise ValueError('--interactive and --pdf are incompatible options')

    # header and data array
    header = reduced_image[0].header
    mecs_header = datamodel.get_mecs_header(reduced_image)
    image2d = reduced_image[0].data

    # check grism and filter
    filter_name = header['filter']
    logger.info('Filter: ' + filter_name)
    if filter_name != bound_param.tags['filter']:
        raise ValueError('Filter name does not match!')
    grism_name = header['grism']
    logger.info('Grism: ' + grism_name)
    if grism_name != bound_param.tags['grism']:
        raise ValueError('Grism name does not match!')

    # read the CSU configuration from the image header
    csu_conf = CsuConfiguration.define_from_header(mecs_header)
    logger.debug(csu_conf)

    # read the DTU configuration from the image header
    dtu2_conf = DtuConf.from_img(reduced_image)
    logger.debug("%s", dtu2_conf)

    # set boundary parameters
    parmodel = bound_param.meta_info['parmodel']
    params = bound_params_from_dict(bound_param.__getstate__())
    if abs(debugplot) >= 10:
        print('-' * 83)
        print('* FITTED BOUND PARAMETERS')
        params.pretty_print()
        pause_debugplot(debugplot)

    # determine parameters according to grism+filter combination
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    islitlet_min = wv_parameters['islitlet_min']
    islitlet_max = wv_parameters['islitlet_max']
    if args_nbrightlines is None:
        nbrightlines = wv_parameters['nbrightlines']
    else:
        nbrightlines = [int(idum) for idum in args_nbrightlines.split(',')]
    poly_crval1_linear = wv_parameters['poly_crval1_linear']
    poly_cdelt1_linear = wv_parameters['poly_cdelt1_linear']
    wvmin_expected = wv_parameters['wvmin_expected']
    wvmax_expected = wv_parameters['wvmax_expected']
    wvmin_useful = wv_parameters['wvmin_useful']
    wvmax_useful = wv_parameters['wvmax_useful']

    # list of slitlets to be computed
    logger.info('list_slitlets: [' + str(islitlet_min) + ',... ' +
                str(islitlet_max) + ']')

    # read master arc line wavelengths (only brightest lines)
    wv_master = read_wv_master_from_array(
        master_table=lines_catalog, lines='brightest', debugplot=debugplot
    )

    # read master arc line wavelengths (whole data set)
    wv_master_all = read_wv_master_from_array(
        master_table=lines_catalog, lines='all', debugplot=debugplot
    )

    # check that the arc lines in the master file are properly sorted
    # in ascending order
    for i in range(len(wv_master_all) - 1):
        if wv_master_all[i] >= wv_master_all[i + 1]:
            logger.error('>>> wavelengths: ' +
                         str(wv_master_all[i]) + '  ' +
                         str(wv_master_all[i+1]))
            raise ValueError('Arc lines are not sorted in master file')

    # ---

    image2d_55sp = np.zeros((EMIR_NBARS, EMIR_NAXIS1))

    # compute rectification transformation and wavelength calibration
    # polynomials

    measured_slitlets = []

    cout = '0'
    for islitlet in range(1, EMIR_NBARS + 1):

        if islitlet_min <= islitlet <= islitlet_max:

            # define Slitlet2dArc object
            slt = Slitlet2dArc(
                islitlet=islitlet,
                csu_conf=csu_conf,
                ymargin_bb=args_ymargin_bb,
                params=params,
                parmodel=parmodel,
                debugplot=debugplot
            )

            # extract 2D image corresponding to the selected slitlet, clipping
            # the image beyond the unrectified slitlet (in order to isolate
            # the arc lines of the current slitlet; otherwise there are
            # problems with arc lines from neighbour slitlets)
            image2d_tmp = select_unrectified_slitlet(
                image2d=image2d,
                islitlet=islitlet,
                csu_bar_slit_center=csu_conf.csu_bar_slit_center(islitlet),
                params=params,
                parmodel=parmodel,
                maskonly=False
            )
            slitlet2d = slt.extract_slitlet2d(image2d_tmp)

            # subtract smooth background computed as follows:
            # - median collapsed spectrum of the whole slitlet2d
            # - independent median filtering of the previous spectrum in the
            #   two halves in the spectral direction
            if args_remove_sp_background:
                spmedian = np.median(slitlet2d, axis=0)
                naxis1_tmp = spmedian.shape[0]
                jmidpoint = naxis1_tmp // 2
                sp1 = medfilt(spmedian[:jmidpoint], [201])
                sp2 = medfilt(spmedian[jmidpoint:], [201])
                spbackground = np.concatenate((sp1, sp2))
                slitlet2d -= spbackground

            # locate unknown arc lines
            slt.locate_unknown_arc_lines(
                slitlet2d=slitlet2d,
                times_sigma_threshold=args_times_sigma_threshold)

            # continue working with current slitlet only if arc lines have
            # been detected
            if slt.list_arc_lines is not None:

                # compute intersections between spectrum trails and arc lines
                slt.xy_spectrail_arc_intersections(slitlet2d=slitlet2d)

                # compute rectification transformation
                slt.estimate_tt_to_rectify(order=args_order_fmap,
                                           slitlet2d=slitlet2d)

                # rectify image
                slitlet2d_rect = slt.rectify(slitlet2d,
                                             resampling=2,
                                             transformation=1)

                # median spectrum and line peaks from rectified image
                sp_median, fxpeaks = slt.median_spectrum_from_rectified_image(
                    slitlet2d_rect,
                    sigma_gaussian_filtering=args_sigma_gaussian_filtering,
                    nwinwidth_initial=5,
                    nwinwidth_refined=5,
                    times_sigma_threshold=5,
                    npix_avoid_border=6,
                    nbrightlines=nbrightlines
                )

                image2d_55sp[islitlet - 1, :] = sp_median

                # determine expected wavelength limits prior to the wavelength
                # calibration
                csu_bar_slit_center = csu_conf.csu_bar_slit_center(islitlet)
                crval1_linear = poly_crval1_linear(csu_bar_slit_center)
                cdelt1_linear = poly_cdelt1_linear(csu_bar_slit_center)
                expected_wvmin = crval1_linear - \
                                 args_margin_npix * cdelt1_linear
                naxis1_linear = sp_median.shape[0]
                crvaln_linear = crval1_linear + \
                                (naxis1_linear - 1) * cdelt1_linear
                expected_wvmax = crvaln_linear + \
                                 args_margin_npix * cdelt1_linear
                # override previous estimates when necessary
                if wvmin_expected is not None:
                    expected_wvmin = wvmin_expected
                if wvmax_expected is not None:
                    expected_wvmax = wvmax_expected

                # clip initial master arc line list with bright lines to
                # the expected wavelength range
                lok1 = expected_wvmin <= wv_master
                lok2 = wv_master <= expected_wvmax
                lok = lok1 * lok2
                wv_master_eff = wv_master[lok]

                # perform initial wavelength calibration
                solution_wv = wvcal_spectrum(
                    sp=sp_median,
                    fxpeaks=fxpeaks,
                    poly_degree_wfit=args_poldeg_initial,
                    wv_master=wv_master_eff,
                    wv_ini_search=expected_wvmin,
                    wv_end_search=expected_wvmax,
                    wvmin_useful=wvmin_useful,
                    wvmax_useful=wvmax_useful,
                    geometry=args_geometry,
                    debugplot=slt.debugplot
                )
                # store initial wavelength calibration polynomial in current
                # slitlet instance
                slt.wpoly = np.polynomial.Polynomial(solution_wv.coeff)
                pause_debugplot(debugplot)

                # clip initial master arc line list with all the lines to
                # the expected wavelength range
                lok1 = expected_wvmin <= wv_master_all
                lok2 = wv_master_all <= expected_wvmax
                lok = lok1 * lok2
                wv_master_all_eff = wv_master_all[lok]

                # clip master arc line list to useful region
                if wvmin_useful is not None:
                    lok = wvmin_useful <= wv_master_all_eff
                    wv_master_all_eff  = wv_master_all_eff[lok]
                if wvmax_useful is not None:
                    lok = wv_master_all_eff <= wvmax_useful
                    wv_master_all_eff  = wv_master_all_eff[lok]

                # refine wavelength calibration
                if args_poldeg_refined > 0:
                    plottitle = '[slitlet#{}, refined]'.format(islitlet)
                    poly_refined, yres_summary = refine_arccalibration(
                        sp=sp_median,
                        poly_initial=slt.wpoly,
                        wv_master=wv_master_all_eff,
                        poldeg=args_poldeg_refined,
                        ntimes_match_wv=1,
                        interactive=args_interactive,
                        threshold=args_threshold_wv,
                        plottitle=plottitle,
                        ylogscale=args_ylogscale,
                        geometry=args_geometry,
                        pdf=args_pdf,
                        debugplot=slt.debugplot
                    )
                    # store refined wavelength calibration polynomial in
                    # current slitlet instance
                    slt.wpoly = poly_refined

                # compute approximate linear values for CRVAL1 and CDELT1
                naxis1_linear = sp_median.shape[0]
                crmin1_linear = slt.wpoly(1)
                crmax1_linear = slt.wpoly(naxis1_linear)
                slt.crval1_linear = crmin1_linear
                slt.cdelt1_linear = \
                    (crmax1_linear - crmin1_linear) / (naxis1_linear - 1)

                # check that the trimming of wv_master and wv_master_all has
                # preserved the wavelength range [crmin1_linear, crmax1_linear]
                if crmin1_linear < expected_wvmin:
                    logger.warning(">>> islitlet: " +str(islitlet))
                    logger.warning("expected_wvmin: " + str(expected_wvmin))
                    logger.warning("crmin1_linear.: " + str(crmin1_linear))
                    logger.warning("WARNING: Unexpected crmin1_linear < "
                                   "expected_wvmin")
                if crmax1_linear > expected_wvmax:
                    logger.warning(">>> islitlet: " +str(islitlet))
                    logger.warning("expected_wvmax: " + str(expected_wvmax))
                    logger.warning("crmax1_linear.: " + str(crmax1_linear))
                    logger.warning("WARNING: Unexpected crmax1_linear > "
                                   "expected_wvmax")

                cout += '.'

            else:

                cout += 'x'

            if islitlet % 10 == 0:
                if cout != 'x':
                    cout = str(islitlet // 10)

            if debugplot != 0:
                pause_debugplot(debugplot)

        else:

            # define Slitlet2dArc object
            slt = Slitlet2dArc(
                islitlet=islitlet,
                csu_conf=csu_conf,
                ymargin_bb=args_ymargin_bb,
                params=None,
                parmodel=None,
                debugplot=debugplot
            )

            cout += 'i'

        # store current slitlet in list of measured slitlets
        measured_slitlets.append(slt)

        logger.info(cout)

    # ---

    # generate FITS file structure with 55 spectra corresponding to the
    # median spectrum for each slitlet
    reduced_55sp_l = copy_img(reduced_image)
    reduced_55sp = reduced_55sp_l[0]
    reduced_55sp.data = image2d_55sp
    reduced_55sp.header['crpix1'] = (0.0, 'reference pixel')
    reduced_55sp.header['crval1'] = (0.0, 'central value at crpix2')
    reduced_55sp.header['cdelt1'] = (1.0, 'increment')
    reduced_55sp.header['ctype1'] = 'PIXEL'
    reduced_55sp.header['cunit1'] = ('Pixel', 'units along axis2')
    reduced_55sp.header['crpix2'] = (0.0, 'reference pixel')
    reduced_55sp.header['crval2'] = (0.0, 'central value at crpix2')
    reduced_55sp.header['cdelt2'] = (1.0, 'increment')
    reduced_55sp.header['ctype2'] = 'PIXEL'
    reduced_55sp.header['cunit2'] = ('Pixel', 'units along axis2')

    # ---

    # Generate structure to store intermediate results
    outdict = {}
    outdict['instrument'] = 'EMIR'
    outdict['meta_info'] = {}
    outdict['meta_info']['creation_date'] = datetime.now().isoformat()
    outdict['meta_info']['description'] = \
        'computation of rectification and wavelength calibration polynomial ' \
        'coefficients for a particular CSU configuration'
    outdict['meta_info']['recipe_name'] = 'undefined'
    outdict['meta_info']['origin'] = {}
    outdict['meta_info']['origin']['bound_param_uuid'] = \
        bound_param.uuid
    outdict['meta_info']['origin']['arc_image_uuid'] = 'undefined'
    outdict['tags'] = {}
    outdict['tags']['grism'] = grism_name
    outdict['tags']['filter'] = filter_name
    outdict['tags']['islitlet_min'] = islitlet_min
    outdict['tags']['islitlet_max'] = islitlet_max
    outdict['dtu_configuration'] = dtu2_conf.outdict()
    outdict['uuid'] = str(uuid4())
    outdict['contents'] = {}

    missing_slitlets = []
    for slt in measured_slitlets:

        islitlet = slt.islitlet

        if islitlet_min <= islitlet <= islitlet_max:

            # avoid error when creating a python list of coefficients from
            # numpy polynomials when the polynomials do not exist (note that
            # the JSON format doesn't handle numpy arrays and such arrays must
            # be transformed into native python lists)
            if slt.wpoly is None:
                wpoly_coeff = None
            else:
                wpoly_coeff = slt.wpoly.coef.tolist()
            if slt.wpoly_longslit_model is None:
                wpoly_coeff_longslit_model = None
            else:
                wpoly_coeff_longslit_model = \
                    slt.wpoly_longslit_model.coef.tolist()

            # avoid similar error when creating a python list of coefficients
            # when the numpy array does not exist; note that this problem
            # does not happen with tt?_aij_longslit_model and
            # tt?_bij_longslit_model because the latter have already been
            # created as native python lists
            if slt.ttd_aij is None:
                ttd_aij = None
            else:
                ttd_aij = slt.ttd_aij.tolist()
            if slt.ttd_bij is None:
                ttd_bij = None
            else:
                ttd_bij = slt.ttd_bij.tolist()
            if slt.tti_aij is None:
                tti_aij = None
            else:
                tti_aij = slt.tti_aij.tolist()
            if slt.tti_bij is None:
                tti_bij = None
            else:
                tti_bij = slt.tti_bij.tolist()

            # creating temporary dictionary with the information corresponding
            # to the current slitlett that will be saved in the JSON file
            tmp_dict = {
                'csu_bar_left': slt.csu_bar_left,
                'csu_bar_right': slt.csu_bar_right,
                'csu_bar_slit_center': slt.csu_bar_slit_center,
                'csu_bar_slit_width': slt.csu_bar_slit_width,
                'x0_reference': slt.x0_reference,
                'y0_reference_lower': slt.y0_reference_lower,
                'y0_reference_middle': slt.y0_reference_middle,
                'y0_reference_upper': slt.y0_reference_upper,
                'y0_reference_lower_expected':
                    slt.y0_reference_lower_expected,
                'y0_reference_middle_expected':
                    slt.y0_reference_middle_expected,
                'y0_reference_upper_expected':
                    slt.y0_reference_upper_expected,
                'y0_frontier_lower': slt.y0_frontier_lower,
                'y0_frontier_upper': slt.y0_frontier_upper,
                'y0_frontier_lower_expected': slt.y0_frontier_lower_expected,
                'y0_frontier_upper_expected': slt.y0_frontier_upper_expected,
                'corr_yrect_a': slt.corr_yrect_a,
                'corr_yrect_b': slt.corr_yrect_b,
                'min_row_rectified': slt.min_row_rectified,
                'max_row_rectified': slt.max_row_rectified,
                'ymargin_bb': slt.ymargin_bb,
                'bb_nc1_orig': slt.bb_nc1_orig,
                'bb_nc2_orig': slt.bb_nc2_orig,
                'bb_ns1_orig': slt.bb_ns1_orig,
                'bb_ns2_orig': slt.bb_ns2_orig,
                'spectrail': {
                    'poly_coef_lower':
                        slt.list_spectrails[
                            slt.i_lower_spectrail].poly_funct.coef.tolist(),
                    'poly_coef_middle':
                        slt.list_spectrails[
                            slt.i_middle_spectrail].poly_funct.coef.tolist(),
                    'poly_coef_upper':
                        slt.list_spectrails[
                            slt.i_upper_spectrail].poly_funct.coef.tolist(),
                },
                'frontier': {
                    'poly_coef_lower':
                        slt.list_frontiers[0].poly_funct.coef.tolist(),
                    'poly_coef_upper':
                        slt.list_frontiers[1].poly_funct.coef.tolist(),
                },
                'ttd_order': slt.ttd_order,
                'ttd_aij': ttd_aij,
                'ttd_bij': ttd_bij,
                'tti_aij': tti_aij,
                'tti_bij': tti_bij,
                'ttd_order_longslit_model': slt.ttd_order_longslit_model,
                'ttd_aij_longslit_model': slt.ttd_aij_longslit_model,
                'ttd_bij_longslit_model': slt.ttd_bij_longslit_model,
                'tti_aij_longslit_model': slt.tti_aij_longslit_model,
                'tti_bij_longslit_model': slt.tti_bij_longslit_model,
                'wpoly_coeff': wpoly_coeff,
                'wpoly_coeff_longslit_model': wpoly_coeff_longslit_model,
                'crval1_linear': slt.crval1_linear,
                'cdelt1_linear': slt.cdelt1_linear
            }
        else:
            missing_slitlets.append(islitlet)
            tmp_dict = {
                'csu_bar_left': slt.csu_bar_left,
                'csu_bar_right': slt.csu_bar_right,
                'csu_bar_slit_center': slt.csu_bar_slit_center,
                'csu_bar_slit_width': slt.csu_bar_slit_width,
                'x0_reference': slt.x0_reference,
                'y0_frontier_lower_expected': slt.y0_frontier_lower_expected,
                'y0_frontier_upper_expected': slt.y0_frontier_upper_expected
            }
        slitlet_label = "slitlet" + str(islitlet).zfill(2)
        outdict['contents'][slitlet_label] = tmp_dict

    # ---

    # OBSOLETE
    '''
    # save JSON file needed to compute the MOS model
    with open(args.out_json.name, 'w') as fstream:
        json.dump(outdict, fstream, indent=2, sort_keys=True)
        print('>>> Saving file ' + args.out_json.name)
    '''

    # ---

    # Create object of type RectWaveCoeff with coefficients for
    # rectification and wavelength calibration
    rectwv_coeff = RectWaveCoeff(instrument='EMIR')
    rectwv_coeff.quality_control = numina.types.qc.QC.GOOD
    rectwv_coeff.tags['grism'] = grism_name
    rectwv_coeff.tags['filter'] = filter_name
    rectwv_coeff.meta_info['origin']['bound_param'] = \
        'uuid' + bound_param.uuid
    rectwv_coeff.meta_info['dtu_configuration'] = outdict['dtu_configuration']
    rectwv_coeff.total_slitlets = EMIR_NBARS
    rectwv_coeff.missing_slitlets = missing_slitlets
    for i in range(EMIR_NBARS):
        islitlet = i + 1
        dumdict = {'islitlet': islitlet}
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        if cslitlet in outdict['contents']:
            dumdict.update(outdict['contents'][cslitlet])
        else:
            raise ValueError("Unexpected error")
        rectwv_coeff.contents.append(dumdict)
    # debugging __getstate__ and __setstate__
    # rectwv_coeff.writeto(args.out_json.name)
    # print('>>> Saving file ' + args.out_json.name)
    # check_setstate_getstate(rectwv_coeff, args.out_json.name)
    logger.info('Generating RectWaveCoeff object with uuid=' +
                rectwv_coeff.uuid)

    return rectwv_coeff, reduced_55sp


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: determine rectification and wavelength '
                    'calibration polynomials from arc image'
    )

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file with longslit data",
                        type=argparse.FileType('rb'))
    parser.add_argument("--bound_param", required=True,
                        help="Input JSON with fitted boundary parameters",
                        type=argparse.FileType('rt'))
    parser.add_argument("--order_fmap", required=True,
                        help="Order of the 2D rectification transformation "
                             "(default=2)",
                        default=2, type=int)
    parser.add_argument("--wv_master_file", required=True,
                        help="TXT file containing wavelengths")
    parser.add_argument("--poldeg_initial", required=True,
                        help="Polynomial degree for initial calibration",
                        type=int)
    parser.add_argument("--poldeg_refined", required=True,
                        help="Polynomial degree for refined calibration "
                             "(0=do not refine)",
                        type=int)
    parser.add_argument("--out_json", required=True,
                        help="Output JSON file with results",
                        type=lambda x: arg_file_is_new(parser, x))

    # optional arguments
    parser.add_argument("--interactive",
                        help="Ask the user for confirmation before updating "
                             "the wavelength calibration polynomial",
                        action="store_true")
    parser.add_argument("--ymargin_bb",
                        help="Number of pixels above and below frontiers to "
                             "determine the vertical bounding box of each "
                             "undistorted slitlet (default=2)",
                        type=int, default=2)
    parser.add_argument("--remove_sp_background",
                        help="Remove background spectrum prior to arc line "
                             "detection",
                        action="store_true")
    parser.add_argument("--times_sigma_threshold",
                        help="Times sigma above threshold to detect unknown"
                             " arc lines (default=10)",
                        type=float, default=10)
    parser.add_argument("--margin_npix",
                        help="Number of pixels before and after expected "
                             "wavelength calibrated spectrum to trim the "
                             "wv_master table in the wavelength direction "
                             "(default=50)",
                        type=int, default=50)
    parser.add_argument("--nbrightlines",
                        help="tuple with number of brightlines to "
                             "be employed in the initial wavelength "
                             "calibration (e.g. \"10,5,4\")")
    parser.add_argument("--threshold_wv",
                        help="Minimum signal in the line peaks (default=0)",
                        default=0, type=float)
    parser.add_argument("--sigma_gaussian_filtering",
                        help="Sigma of the gaussian filter to be applied to "
                             "the spectrum in order to avoid problems with "
                             "saturated lines in the wavelength calibration "
                             "process",
                        default=0, type=float)
    parser.add_argument("--out_55sp",
                        help="FITS file containing the set of averaged "
                             "spectra employed to derive the wavelength "
                             "calibration",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))
    parser.add_argument("--ylogscale",
                        help="Display spectrum signal in logarithmic units",
                        action="store_true")
    parser.add_argument("--geometry",
                        help="tuple x,y,dx,dy (default 0,0,640,480)",
                        default="0,0,640,480")
    parser.add_argument("--pdffile",
                        help="output PDF file name",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting & debugging options"
                             " (default=0)",
                        default=0, type=int,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")
    args = parser.parse_args(args)

    if args.echo:
        print('\033[1m\033[31m% ' + ' '.join(sys.argv) + '\033[0m\n')

    # ---

    logger = logging.getLogger(__name__)

    logging_from_debugplot(args.debugplot)

    # read pdffile
    if args.pdffile is not None:
        if args.interactive:
            raise ValueError('--interactive is not compatible with --pdffile')
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(args.pdffile.name)
    else:
        pdf = None

    # geometry
    if args.geometry is None:
        geometry = None
    else:
        tmp_str = args.geometry.split(",")
        x_geom = int(tmp_str[0])
        y_geom = int(tmp_str[1])
        dx_geom = int(tmp_str[2])
        dy_geom = int(tmp_str[3])
        geometry = x_geom, y_geom, dx_geom, dy_geom

    # generate HDUList object
    hdulist = fits.open(args.fitsfile)

    # generate RefinedBoundaryModelParam object
    bound_param = RefinedBoundaryModelParam._datatype_load(
        args.bound_param.name)

    # generate lines_catalog
    lines_catalog = np.genfromtxt(args.wv_master_file)

    rectwv_coeff, reduced_55sp = rectwv_coeff_from_arc_image(
        hdulist,
        bound_param,
        lines_catalog,
        args_nbrightlines=args.nbrightlines,
        args_ymargin_bb=args.ymargin_bb,
        args_remove_sp_background=args.remove_sp_background,
        args_times_sigma_threshold=args.times_sigma_threshold,
        args_order_fmap=args.order_fmap,
        args_sigma_gaussian_filtering=args.sigma_gaussian_filtering,
        args_margin_npix=args.margin_npix,
        args_poldeg_initial=args.poldeg_initial,
        args_poldeg_refined=args.poldeg_refined,
        args_interactive=args.interactive,
        args_threshold_wv=args.threshold_wv,
        args_ylogscale=args.ylogscale,
        args_pdf=pdf,
        args_geometry=geometry,
        debugplot=args.debugplot
    )

    # save image with collapsed spectra employed to determine the
    # wavelength calibration
    if args.out_55sp is not None:
        reduced_55sp.writeto(args.out_55sp, overwrite=True)

    # save RectWaveCoeff object into JSON file
    rectwv_coeff.writeto(args.out_json.name)
    logger.info('>>> Saving file ' + args.out_json.name)
    # debugging __getstate__ and __setstate__
    # check_setstate_getstate(rectwv_coeff, args.out_json.name)

    if pdf is not None:
        pdf.close()


if __name__ == "__main__":
    main()