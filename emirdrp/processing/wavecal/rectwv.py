#
# Copyright 2018 Universidad Complutense de Madrid
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


"""Wavelength calibration functionality"""

from __future__ import division, print_function

from astropy.io import fits
from datetime import datetime
import logging
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import medfilt
import sys
import time
from uuid import uuid4

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.distortion import ncoef_fmap
from numina.array.wavecalib.__main__ import read_wv_master_from_array
from numina.array.wavecalib.__main__ import wvcal_spectrum
from numina.array.wavecalib.arccalibration import refine_arccalibration
from numina.array.wavecalib.resample import resample_image2d_flux
import numina.types.qc

from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.instrument.dtu_configuration import DtuConfiguration
from emirdrp.products import RectWaveCoeff
from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.fit_boundaries import expected_distorted_boundaries
from emirdrp.tools.fit_boundaries import expected_distorted_frontiers
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers
from emirdrp.tools.select_unrectified_slitlets import \
    select_unrectified_slitlet

from .islitlet_progress import islitlet_progress
from .set_wv_parameters import set_wv_parameters
from .slitlet2darc import Slitlet2dArc
from .slitlet2d import Slitlet2D

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
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

    # header and data array
    header = reduced_image[0].header
    image2d = reduced_image[0].data

    # check grism and filter
    filter_name = header['filter']
    logger.debug('Filter: ' + filter_name)
    if filter_name != bound_param.tags['filter']:
        raise ValueError('Filter name does not match!')
    grism_name = header['grism']
    logger.debug('Grism: ' + grism_name)
    if grism_name != bound_param.tags['grism']:
        raise ValueError('Grism name does not match!')

    # read the CSU configuration from the image header
    csu_conf = CsuConfiguration.define_from_header(header)
    logger.debug(csu_conf)

    # read the DTU configuration from the image header
    dtu_conf = DtuConfiguration.define_from_header(header)
    logger.debug(dtu_conf)

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

    # list of slitlets to be computed
    list_slitlets = list(range(islitlet_min, islitlet_max + 1))
    print('list_slitlets:\n', list_slitlets)

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
            print('>>> wavelengths: ', wv_master_all[i], wv_master_all[i + 1])
            raise ValueError('Arc lines are not sorted in master file')

    # ---

    image2d_55sp = np.zeros((EMIR_NBARS, EMIR_NAXIS1))

    # compute rectification transformation and wavelength calibration
    # polynomials

    measured_slitlets = []

    for islitlet in list_slitlets:

        # define Slitlet2dArc object
        slt = Slitlet2dArc(
            islitlet=islitlet,
            params=params,
            parmodel=parmodel,
            csu_conf=csu_conf,
            ymargin_bb=args_ymargin_bb,
            debugplot=debugplot
        )

        # extract 2D image corresponding to the selected slitlet, clipping
        # the image beyond the unrectified slitlet (in order to isolate
        # the arc lines of the current slitlet; otherwise there are problems
        # with arc lines from neighbour slitlets)
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
            expected_wvmin = crval1_linear - args_margin_npix * cdelt1_linear
            naxis1_linear = sp_median.shape[0]
            crvaln_linear = crval1_linear + (naxis1_linear - 1) * cdelt1_linear
            expected_wvmax = crvaln_linear + args_margin_npix * cdelt1_linear

            # clip initial master arc line list with bright lines to expected
            # wavelength range
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
                geometry=args_geometry,
                debugplot=slt.debugplot
            )
            # store initial wavelength calibration polynomial in current
            # slitlet instance
            slt.wpoly = np.polynomial.Polynomial(solution_wv.coeff)
            pause_debugplot(debugplot)

            # clip initial master arc line list with all the lines to expected
            # wavelength range
            lok1 = expected_wvmin <= wv_master_all
            lok2 = wv_master_all <= expected_wvmax
            lok = lok1 * lok2
            wv_master_all_eff = wv_master_all[lok]

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
                # store refined wavelength calibration polynomial in current
                # slitlet instance
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
                print(">>> islitlet: ", islitlet)
                print(">>> expected_wvmin:", expected_wvmin)
                print(">>> crmin1_linear.:", crmin1_linear)
                print(">>> WARNING: Unexpected crmin1_linear < expected_wvmin")
            if crmax1_linear > expected_wvmax:
                print(">>> islitlet: ", islitlet)
                print(">>> expected_wvmax:", expected_wvmax)
                print(">>> crmax1_linear.:", crmax1_linear)
                print(">>> WARNING: Unexpected crmax1_linear > expected_wvmin")

            cout = '.'

        else:

            cout = 'x'

        # store current slitlet in list of measured slitlets
        measured_slitlets.append(slt)

        if debugplot == 0:
            if islitlet % 10 == 0:
                if cout != 'x':
                    cout = str(islitlet // 10)
            sys.stdout.write(cout)
            # print(slt)
            if islitlet == list_slitlets[-1]:
                sys.stdout.write('\n')
            sys.stdout.flush()
        else:
            pause_debugplot(debugplot)

    # ---

    # ToDo: return out_55sp ???
    '''
    # Save image with collapsed spectra employed to determine the
    # wavelength calibration
    if args.out_55sp is not None:
        save_ndarray_to_fits(
            array=image2d_55sp,
            file_name=args.out_55sp,
            cast_to_float=True,
            overwrite=True
        )
    '''
    reduced_55sp = fits.PrimaryHDU(data=image2d_55sp)
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
    outdict['meta-info'] = {}
    outdict['meta-info']['creation_date'] = datetime.now().isoformat()
    outdict['meta-info']['description'] = \
        'computation of rectification and wavelength calibration polynomial ' \
        'coefficients for a particular CSU configuration'
    outdict['meta-info']['recipe_name'] = 'undefined'
    outdict['meta-info']['origin'] = {}
    outdict['meta-info']['origin']['bound_param_uuid'] = \
        bound_param.uuid
    outdict['meta-info']['origin']['arc_image_uuid'] = 'undefined'
    outdict['tags'] = {}
    outdict['tags']['grism'] = grism_name
    outdict['tags']['filter'] = filter_name
    outdict['tags']['islitlet_min'] = islitlet_min
    outdict['tags']['islitlet_max'] = islitlet_max
    outdict['dtu_configuration'] = dtu_conf.outdict()
    outdict['uuid'] = str(uuid4())
    outdict['contents'] = {}

    for slt in measured_slitlets:

        islitlet = slt.islitlet

        # avoid error when creating a python list of coefficients from
        # numpy polynomials when the polynomials do not exist (note that the
        # JSON format doesn't handle numpy arrays and such arrays must be
        # transformed into native python lists)
        if slt.wpoly is None:
            wpoly_coeff = None
        else:
            wpoly_coeff = slt.wpoly.coef.tolist()
        if slt.wpoly_longslit_model is None:
            wpoly_coeff_longslit_model = None
        else:
            wpoly_coeff_longslit_model = slt.wpoly_longslit_model.coef.tolist()

        # avoid similar error when creating a python list of coefficients
        # when the numpy array does not exist; note that this problem
        # does not happen with tt?_aij_longslit_model and
        # tt?_bij_longslit_model because the latter have already been created
        # as native python lists
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

        # creating temporary dictionary with the information corresponding to
        # the current slitlett that will be saved in the JSON file
        tmp_dict = {
            'csu_bar_left': slt.csu_bar_left,
            'csu_bar_right': slt.csu_bar_right,
            'csu_bar_slit_center': slt.csu_bar_slit_center,
            'csu_bar_slit_width': slt.csu_bar_slit_width,
            'x0_reference': slt.x0_reference,
            'y0_reference_lower': slt.y0_reference_lower,
            'y0_reference_middle': slt.y0_reference_middle,
            'y0_reference_upper': slt.y0_reference_upper,
            'y0_frontier_lower': slt.y0_frontier_lower,
            'y0_frontier_upper': slt.y0_frontier_upper,
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
    for i in range(EMIR_NBARS):
        islitlet = i + 1
        dumdict = {'islitlet': islitlet}
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        if cslitlet in outdict['contents']:
            dumdict.update(outdict['contents'][cslitlet])
        else:
            dumdict.update({
                'ttd_order': 0
            })
            rectwv_coeff.missing_slitlets.append(islitlet)
        rectwv_coeff.contents.append(dumdict)
    # debugging __getstate__ and __setstate__
    # rectwv_coeff.writeto(args.out_json.name)
    # print('>>> Saving file ' + args.out_json.name)
    # check_setstate_getstate(rectwv_coeff, args.out_json.name)
    logger.info('Generating RectWaveCoeff object with uuid=' +
                rectwv_coeff.uuid)

    return rectwv_coeff, reduced_55sp


def rectwv_coeff_from_mos_library(reduced_image, master_rectwv,
                                           ignore_DTUconf=True):
    """Evaluate rect.+wavecal. coefficients from MOS library

    Parameters
    ----------
    reduced_image : HDUList object
        Image with preliminary basic reduction: bpm, bias, dark and
        flatfield.
    master_rectwv : MasterRectWave instance
        Rectification and Wavelength Calibrartion Library product.
        Contains the library of polynomial coefficients necessary
        to generate an instance of RectWaveCoeff with the rectification
        and wavelength calibration coefficients for the particular
        CSU configuration.
    ignore_DTUconf : bool
        If True, ignore differences in DTU configuration.

    Returns
    -------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.

    """

    logger = logging.getLogger(__name__)

    # header
    header = reduced_image[0].header

    # read the CSU configuration from the image header
    csu_conf = CsuConfiguration.define_from_header(header)
    logger.debug(csu_conf)

    # read the DTU configuration from the image header
    dtu_conf = DtuConfiguration.define_from_header(header)
    logger.debug(dtu_conf)

    # retrieve DTU configuration from MasterRectWave object
    dtu_conf_calib = DtuConfiguration.define_from_dictionary(
        master_rectwv.meta_info['dtu_configuration']
    )
    # check that the DTU configuration employed to obtain the calibration
    # corresponds to the DTU configuration in the input FITS file
    if dtu_conf != dtu_conf_calib:
        logger.info('DTU configuration from image header:')
        logger.info(dtu_conf)
        logger.info('DTU configuration from master calibration:')
        logger.info(dtu_conf_calib)
        if ignore_DTUconf:
            logger.warning('DTU configuration differences found!')
        else:
            raise ValueError("DTU configurations do not match!")
    else:
        logger.info('DTU configuration match!')

    # check grism and filter
    filter_name = header['filter']
    logger.debug('Filter: ' + filter_name)
    if filter_name != master_rectwv.tags['filter']:
        raise ValueError('Filter name does not match!')
    grism_name = header['grism']
    logger.debug('Grism: ' + grism_name)
    if grism_name != master_rectwv.tags['grism']:
        raise ValueError('Grism name does not match!')

    # valid slitlet numbers
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in master_rectwv.missing_slitlets:
        list_valid_islitlets.remove(idel)
    logger.debug('valid slitlet numbers: ' + str(list_valid_islitlets))

    # initialize intermediate dictionary with relevant information
    # (note: this dictionary corresponds to an old structure employed to
    # store the information in a JSON file; this is no longer necessary,
    # but here we reuse that dictionary for convenience)
    outdict = {}
    outdict['instrument'] = 'EMIR'
    outdict['meta-info'] = {}
    outdict['meta-info']['creation_date'] = datetime.now().isoformat()
    outdict['meta-info']['description'] = \
        'computation of rectification and wavelength calibration polynomial ' \
        'coefficients for a particular CSU configuration from a MOS model '
    outdict['meta-info']['recipe_name'] = 'undefined'
    outdict['meta-info']['origin'] = {}
    outdict['meta-info']['origin']['fits_frame_uuid'] = 'TBD'
    outdict['meta-info']['origin']['rect_wpoly_mos_uuid'] = \
        master_rectwv.uuid
    outdict['meta-info']['origin']['fitted_boundary_param_uuid'] = \
        master_rectwv.meta_info['origin']['bound_param']
    outdict['tags'] = {}
    outdict['tags']['grism'] = grism_name
    outdict['tags']['filter'] = filter_name
    outdict['dtu_configuration'] = dtu_conf.outdict()
    outdict['uuid'] = str(uuid4())
    outdict['contents'] = {}

    # compute rectification and wavelength calibration coefficients for each
    # slitlet according to its csu_bar_slit_center value
    for islitlet in list_valid_islitlets:
        cslitlet = 'slitlet' + str(islitlet).zfill(2)

        # csu_bar_slit_center of current slitlet in initial FITS image
        csu_bar_slit_center = csu_conf.csu_bar_slit_center(islitlet)

        # input data structure
        tmpdict = master_rectwv.contents[islitlet - 1]
        list_csu_bar_slit_center = tmpdict['list_csu_bar_slit_center']

        # check extrapolations
        if csu_bar_slit_center < min(list_csu_bar_slit_center):
            logger.warning('extrapolating table with ' + cslitlet)
            logger.warning('minimum tabulated value: ' +
                           str(min(list_csu_bar_slit_center)))
            logger.warning('sought value...........: ' +
                           str(csu_bar_slit_center))
        if csu_bar_slit_center > max(list_csu_bar_slit_center):
            logger.warning('extrapolating table with ' + cslitlet)
            logger.warning('maximum tabulated value: ' +
                           str(max(list_csu_bar_slit_center)))
            logger.warning('sought value...........: ' +
                           str(csu_bar_slit_center))

        # rectification coefficients
        ttd_order = tmpdict['ttd_order']
        ncoef = ncoef_fmap(ttd_order)
        outdict['contents'][cslitlet] = {}
        outdict['contents'][cslitlet]['ttd_order'] = ttd_order
        outdict['contents'][cslitlet]['ttd_order_longslit_model'] = None
        for keycoef in ['ttd_aij', 'ttd_bij', 'tti_aij', 'tti_bij']:
            coef_out = []
            for icoef in range(ncoef):
                ccoef = str(icoef).zfill(2)
                list_cij = tmpdict['list_' + keycoef + '_' + ccoef]
                funinterp_coef = interp1d(list_csu_bar_slit_center,
                                          list_cij,
                                          kind='linear',
                                          fill_value='extrapolate')
                # note: funinterp_coef expects a numpy array
                dum = funinterp_coef([csu_bar_slit_center])
                coef_out.append(dum[0])
            outdict['contents'][cslitlet][keycoef] = coef_out
            outdict['contents'][cslitlet][keycoef + '_longslit_model'] = None

        # wavelength calibration coefficients
        ncoef = tmpdict['wpoly_degree'] + 1
        wpoly_coeff = []
        for icoef in range(ncoef):
            ccoef = str(icoef).zfill(2)
            list_cij = tmpdict['list_wpoly_coeff_' + ccoef]
            funinterp_coef = interp1d(list_csu_bar_slit_center,
                                      list_cij,
                                      kind='linear',
                                      fill_value='extrapolate')
            # note: funinterp_coef expects a numpy array
            dum = funinterp_coef([csu_bar_slit_center])
            wpoly_coeff.append(dum[0])
        outdict['contents'][cslitlet]['wpoly_coeff'] = wpoly_coeff
        outdict['contents'][cslitlet]['wpoly_coeff_longslit_model'] = None

        # update cdelt1_linear and crval1_linear
        wpoly_function = np.polynomial.Polynomial(wpoly_coeff)
        crmin1_linear = wpoly_function(1)
        crmax1_linear = wpoly_function(EMIR_NAXIS1)
        cdelt1_linear = (crmax1_linear - crmin1_linear) / (EMIR_NAXIS1 - 1)
        crval1_linear = crmin1_linear
        outdict['contents'][cslitlet]['crval1_linear'] = crval1_linear
        outdict['contents'][cslitlet]['cdelt1_linear'] = cdelt1_linear

        # update CSU keywords
        outdict['contents'][cslitlet]['csu_bar_left'] = \
            csu_conf.csu_bar_left(islitlet)
        outdict['contents'][cslitlet]['csu_bar_right'] = \
            csu_conf.csu_bar_right(islitlet)
        outdict['contents'][cslitlet]['csu_bar_slit_center'] = \
            csu_conf.csu_bar_slit_center(islitlet)
        outdict['contents'][cslitlet]['csu_bar_slit_width'] = \
            csu_conf.csu_bar_slit_width(islitlet)

    # for each slitlet compute spectrum trails and frontiers using the
    # fitted boundary parameters
    fitted_bound_param_json = {
        'contents': master_rectwv.meta_info['refined_boundary_model']
    }
    parmodel = fitted_bound_param_json['contents']['parmodel']
    fitted_bound_param_json.update({'meta_info': {'parmodel': parmodel}})
    params = bound_params_from_dict(fitted_bound_param_json)
    logger.debug('Fitted boundary parameters:')
    logger.debug(params.pretty_print())
    for islitlet in list_valid_islitlets:
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        # csu_bar_slit_center of current slitlet in initial FITS image
        csu_bar_slit_center = csu_conf.csu_bar_slit_center(islitlet)
        # compute and store x0_reference value
        x0_reference = float(EMIR_NAXIS1) / 2.0 + 0.5
        outdict['contents'][cslitlet]['x0_reference'] = x0_reference
        # compute spectrum trails (lower, middle and upper)
        list_spectrails = expected_distorted_boundaries(
            islitlet, csu_bar_slit_center,
            [0, 0.5, 1], params, parmodel,
            numpts=101, deg=5, debugplot=0
        )
        # store spectrails in output JSON file
        outdict['contents'][cslitlet]['spectrail'] = {}
        for idum, cdum in zip(range(3), ['lower', 'middle', 'upper']):
            outdict['contents'][cslitlet]['spectrail']['poly_coef_' + cdum] = \
                list_spectrails[idum].poly_funct.coef.tolist()
            outdict['contents'][cslitlet]['y0_reference_' + cdum] = \
                list_spectrails[idum].poly_funct(x0_reference)
        # compute frontiers (lower, upper)
        list_frontiers = expected_distorted_frontiers(
            islitlet, csu_bar_slit_center,
            params, parmodel,
            numpts=101, deg=5, debugplot=0
        )
        # store frontiers in output JSON
        outdict['contents'][cslitlet]['frontier'] = {}
        for idum, cdum in zip(range(2), ['lower', 'upper']):
            outdict['contents'][cslitlet]['frontier']['poly_coef_' + cdum] = \
                list_frontiers[idum].poly_funct.coef.tolist()
            outdict['contents'][cslitlet]['y0_frontier_' + cdum] = \
                list_frontiers[idum].poly_funct(x0_reference)

    # store bounding box parameters for each slitlet
    xdum = np.linspace(1, EMIR_NAXIS1, num=EMIR_NAXIS1)
    for islitlet in list_valid_islitlets:
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        # parameters already available in the input JSON file
        for par in ['bb_nc1_orig', 'bb_nc2_orig', 'ymargin_bb']:
            outdict['contents'][cslitlet][par] = \
                master_rectwv.contents[islitlet - 1][par]
        # estimate bb_ns1_orig and bb_ns2_orig using the already computed
        # frontiers and the value of ymargin_bb, following the same approach
        # employed in the script rect_wpoly_from_longslit; see
        # Slitlet2dLongSlitArc.__init__()
        poly_lower_frontier = np.polynomial.Polynomial(
            outdict['contents'][cslitlet]['frontier']['poly_coef_lower']
        )
        poly_upper_frontier = np.polynomial.Polynomial(
            outdict['contents'][cslitlet]['frontier']['poly_coef_upper']
        )
        ylower = poly_lower_frontier(xdum)
        yupper = poly_upper_frontier(xdum)
        ymargin_bb = master_rectwv.contents[islitlet - 1]['ymargin_bb']
        bb_ns1_orig = int(ylower.min() + 0.5) - ymargin_bb
        if bb_ns1_orig < 1:
            bb_ns1_orig = 1
        bb_ns2_orig = int(yupper.max() + 0.5) + ymargin_bb
        if bb_ns2_orig > EMIR_NAXIS2:
            bb_ns2_orig = EMIR_NAXIS2
        outdict['contents'][cslitlet]['bb_ns1_orig'] = bb_ns1_orig
        outdict['contents'][cslitlet]['bb_ns2_orig'] = bb_ns2_orig

    # ---

    # Create object of type RectWaveCoeff with coefficients for
    # rectification and wavelength calibration
    rectwv_coeff = RectWaveCoeff(instrument='EMIR')
    rectwv_coeff.quality_control = numina.types.qc.QC.GOOD
    rectwv_coeff.tags['grism'] = grism_name
    rectwv_coeff.tags['filter'] = filter_name
    rectwv_coeff.meta_info['origin']['bound_param'] = \
        master_rectwv.meta_info['origin']['bound_param']
    rectwv_coeff.meta_info['origin']['master_rectwv'] = \
        'uuid' + master_rectwv.uuid
    rectwv_coeff.meta_info['dtu_configuration'] = outdict['dtu_configuration']
    rectwv_coeff.total_slitlets = EMIR_NBARS
    for i in range(EMIR_NBARS):
        islitlet = i + 1
        dumdict = {'islitlet': islitlet}
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        if cslitlet in outdict['contents']:
            dumdict.update(outdict['contents'][cslitlet])
        else:
            dumdict.update({
                'ttd_order': 0
            })
            rectwv_coeff.missing_slitlets.append(islitlet)
        rectwv_coeff.contents.append(dumdict)
    # debugging __getstate__ and __setstate__
    # rectwv_coeff.writeto(args.out_rect_wpoly.name)
    # print('>>> Saving file ' + args.out_rect_wpoly.name)
    # check_setstate_getstate(rectwv_coeff, args.out_rect_wpoly.name)
    logger.info('Generating RectWaveCoeff object with uuid=' +
                rectwv_coeff.uuid)

    return rectwv_coeff


def apply_rectwv_coeff(reduced_image, rectwv_coeff,
                       resampling=2, ignore_DTUconf=True,
                       debugplot=0):
    """Compute rectification and wavelength calibration coefficients.

        Parameters
        ----------
        reduced_image : HDUList object
            Image with preliminary basic reduction: bpm, bias, dark and
            flatfield.
        rectwv_coeff : RectWaveCoeff instance
            Rectification and wavelength calibration coefficients for the
            particular CSU configuration.
        resampling : int
            1: nearest neighbour, 2: flux preserving interpolation.
        ignore_DTUconf : bool
            If True, ignore differences in DTU configuration.
        debugplot : int
            Debugging level for messages and plots. For details see
            'numina.array.display.pause_debugplot.py'.

        Returns
        -------
        rectwv_image : HDUList object
            Rectified and wavelength calibrated image.

        """

    logger = logging.getLogger(__name__)

    # header and data array
    header = reduced_image[0].header
    image2d = reduced_image[0].data

    # check grism and filter
    filter_name = header['filter']
    logger.debug('Filter: ' + filter_name)
    if filter_name != rectwv_coeff.tags['filter']:
        raise ValueError('Filter name does not match!')
    grism_name = header['grism']
    logger.debug('Grism: ' + grism_name)
    if grism_name != rectwv_coeff.tags['grism']:
        raise ValueError('Grism name does not match!')

    # read the DTU configuration from the image header
    dtu_conf = DtuConfiguration.define_from_header(header)
    logger.debug(dtu_conf)

    # retrieve DTU configuration from RectWaveCoeff object
    dtu_conf_calib = DtuConfiguration.define_from_dictionary(
        rectwv_coeff.meta_info['dtu_configuration']
    )
    # check that the DTU configuration employed to obtain the calibration
    # corresponds to the DTU configuration in the input FITS file
    if dtu_conf != dtu_conf_calib:
        logger.info('DTU configuration from image header:')
        logger.info(dtu_conf)
        logger.info('DTU configuration from master calibration:')
        logger.info(dtu_conf_calib)
        if ignore_DTUconf:
            logger.warning('DTU configuration differences found!')
        else:
            raise ValueError("DTU configurations do not match!")
    else:
        logger.info('DTU configuration match!')

    # valid slitlet numbers
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in rectwv_coeff.missing_slitlets:
        list_valid_islitlets.remove(idel)
    logger.debug('Valid slitlet numbers: ' + str(list_valid_islitlets))

    # ---

    # relevant wavelength calibration parameters for rectified and wavelength
    # calibrated image
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    crpix1_enlarged = wv_parameters['crpix1_enlarged']
    crval1_enlarged = wv_parameters['crval1_enlarged']
    cdelt1_enlarged = wv_parameters['cdelt1_enlarged']
    naxis1_enlarged = wv_parameters['naxis1_enlarged']

    # initialize rectified and wavelength calibrated image
    image2d_rectwv = np.zeros((EMIR_NAXIS2, naxis1_enlarged))

    # main loop
    logger.info('Computing rectification and wavelength calibration')
    for islitlet in list_valid_islitlets:
        if abs(debugplot) >= 10:
            if islitlet == list_valid_islitlets[0]:
                print(time.ctime())
            islitlet_progress(islitlet, EMIR_NBARS)
            if islitlet == list_valid_islitlets[-1]:
                print(' ')
                print(time.ctime())

        # define Slitlet2D object
        slt = Slitlet2D(islitlet=islitlet,
                        rectwv_coeff=rectwv_coeff,
                        debugplot=debugplot)

        # extract (distorted) slitlet from the initial image
        slitlet2d = slt.extract_slitlet2d(image2d)

        # rectify slitlet
        slitlet2d_rect = slt.rectify(slitlet2d, resampling=resampling)

        # wavelength calibration of the rectifed slitlet
        slitlet2d_rect_wv = resample_image2d_flux(
            image2d_orig=slitlet2d_rect,
            naxis1=naxis1_enlarged,
            cdelt1=cdelt1_enlarged,
            crval1=crval1_enlarged,
            crpix1=crpix1_enlarged,
            coeff=slt.wpoly
        )

        # minimum and maximum useful scan (pixel in the spatial direction)
        # for the rectified slitlet
        nscan_min, nscan_max = nscan_minmax_frontiers(
            slt.y0_frontier_lower,
            slt.y0_frontier_upper,
            resize=False
        )
        ii1 = nscan_min - slt.bb_ns1_orig
        ii2 = nscan_max - slt.bb_ns1_orig + 1
        i1 = slt.bb_ns1_orig - 1 + ii1
        i2 = i1 + ii2 - ii1
        image2d_rectwv[i1:i2, :] = slitlet2d_rect_wv[ii1:ii2, :]

        # include scan range in FITS header
        header['sltmin' + str(islitlet).zfill(2)] = i1
        header['sltmax' + str(islitlet).zfill(2)] = i2 - 1

    # modify upper limit of previous slitlet in case of overlapping:
    # note that the overlapped scans have been overwritten with the
    # information from the current slitlet!
    for islitlet in list_valid_islitlets:
        cprevious = 'SLTMAX' + str(islitlet - 1).zfill(2)
        if cprevious in header.keys():
            sltmax_previous = header[cprevious]
            cslitlet = 'SLTMIN' + str(islitlet).zfill(2)
            sltmin_current = header[cslitlet]
            if sltmax_previous >= sltmin_current:
                print('WARNING: ' + cslitlet + '=' +
                      str(sltmin_current).zfill(4) +
                      ' overlaps with ' + cprevious + '=' +
                      str(sltmax_previous).zfill(4) + ' ==> ' + cslitlet +
                      ' set to ' + str(sltmin_current - 1).zfill(4))
                header[cprevious] = sltmin_current - 1

    logger.info('Updating image header')
    # update wavelength calibration in FITS header
    header.remove('crval1')
    header.remove('crpix1')
    header.remove('crval2')
    header.remove('crpix2')
    header['crpix1'] = (crpix1_enlarged, 'reference pixel')
    header['crval1'] = (crval1_enlarged, 'central wavelength at crpix1')
    header['cdelt1'] = (cdelt1_enlarged, 'linear dispersion (Angstrom/pixel)')
    header['cunit1'] = ('Angstrom', 'units along axis1')
    header['ctype1'] = 'WAVELENGTH'
    header['crpix2'] = (0.0, 'reference pixel')
    header['crval2'] = (0.0, 'central value at crpix2')
    header['cdelt2'] = (1.0, 'increment')
    header['ctype2'] = 'PIXEL'
    header['cunit2'] = ('Pixel', 'units along axis2')
    header.remove('cd1_1')
    header.remove('cd1_2')
    header.remove('cd2_1')
    header.remove('cd2_2')
    header.remove('PCD1_1')
    header.remove('PCD1_2')
    header.remove('PCD2_1')
    header.remove('PCD2_2')
    header.remove('PCRPIX1')
    header.remove('PCRPIX2')
    header['history'] = 'Boundary parameters uuid:' + \
                        rectwv_coeff.meta_info['origin']['bound_param'][4:]
    if 'master_rectwv' in rectwv_coeff.meta_info['origin']:
        header['history'] = \
            'MasterRectWave uuid:' + \
            rectwv_coeff.meta_info['origin']['master_rectwv'][4:]
    header['history'] = 'RectWaveCoeff uuid:' + rectwv_coeff.uuid
    header['history'] = 'Rectification and wavelength calibration time ' \
                        + datetime.now().isoformat()

    logger.info('Generating rectified and wavelength calibrated image')

    rectwv_image = fits.PrimaryHDU(data=image2d_rectwv, header=header)

    return rectwv_image
