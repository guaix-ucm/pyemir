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

import argparse
from astropy.io import fits
from datetime import datetime
import numpy as np
from scipy.signal import medfilt
import sys
from uuid import uuid4

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow
from numina.array.wavecalib.__main__ import read_wv_master_file
from numina.array.wavecalib.__main__ import wvcal_spectrum
from numina.array.wavecalib.arccalibration import refine_arccalibration
import numina.types.qc
from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.instrument.dtu_configuration import DtuConfiguration
from emirdrp.processing.wavecal import set_wv_parameters
from emirdrp.processing.wavecal import Slitlet2dArc
from emirdrp.products import RefinedBoundaryModelParam
from emirdrp.products import RectWaveCoeff
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers
from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.save_ndarray_to_fits import save_ndarray_to_fits
from emirdrp.tools.select_unrectified_slitlets import select_unrectified_slitlet

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from numina.tools.arg_file_is_new import arg_file_is_new
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_VALID_FILTERS
from emirdrp.core import EMIR_VALID_GRISMS


def interpolate_bad_rows(image2d):
    """Interpolate quadrant change in original frames.

    Parameters
    ----------
    image2d : numpy array
        Initial image.

    Returns
    -------
    image2d_interpolated : numpy array
        Interpolated image.

    """

    # ToDo: these numbers depend on EMIR_NAXIS1 and EMIR_NAXIS2
    image2d_interpolated = np.copy(image2d)
    image2d_interpolated[1024, :1024] = (image2d[1023, :1024] +
                                         image2d[1025, :1024]) / 2
    image2d_interpolated[1023, 1024:] = (image2d[1022, 1024:] +
                                         image2d[1024, 1024:]) / 2

    return image2d_interpolated


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file with longslit data",
                        type=argparse.FileType('rb'))
    parser.add_argument("--fitted_bound_param", required=True,
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
    parser.add_argument("--out_rect",
                        help="Rectified but not wavelength calibrated output "
                             "FITS file",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))
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

    # read the CSU configuration from the initial FITS file
    csu_conf = CsuConfiguration.define_from_fits(args.fitsfile)
    if abs(args.debugplot) >= 10:
        print(csu_conf)
        pause_debugplot(args.debugplot)

    # read the DTU configuration from the initial FITS file
    dtu_conf = DtuConfiguration.define_from_fits(args.fitsfile)
    if abs(args.debugplot) >= 10:
        print(dtu_conf)
        pause_debugplot(args.debugplot)

    # read fitted boundary parameters
    bound_param = RefinedBoundaryModelParam._datatype_load(
        args.bound_param.name)
    parmodel = bound_param.meta_info['parmodel']
    params = bound_params_from_dict(bound_param.__getstate__())
    if abs(args.debugplot) >= 10:
        print('-' * 83)
        print('* FITTED BOUND PARAMETERS')
        params.pretty_print()
        pause_debugplot(args.debugplot)

    # read input FITS image
    with fits.open(args.fitsfile) as hdulist:
        image2d_header = hdulist[0].header
        image2d = interpolate_bad_rows(hdulist[0].data)

    # read filter and grism names
    grism_name = bound_param.tags['grism']
    filter_name = bound_param.tags['filter']
    if filter_name not in EMIR_VALID_FILTERS:
        raise ValueError('Unexpected filter_name:', filter_name)
    if grism_name not in EMIR_VALID_GRISMS:
        raise ValueError('Unexpected grism_name:', grism_name)
    if filter_name != image2d_header['filter']:
        print('Filter name (FITS file)..:', image2d_header['filter'])
        print('Filter name (bound_param):', filter_name)
        raise ValueError('Filter_name does not match')
    if grism_name != image2d_header['grism']:
        print('Grism name (FITS file)..:', image2d_header['grism'])
        print('Grism name (bound_param):', grism_name)
        raise ValueError('Grism_name does not match')

    # determine parameters according to grism+filter combination
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    islitlet_min = wv_parameters['islitlet_min']
    islitlet_max = wv_parameters['islitlet_max']
    if args.nbrightlines is None:
        nbrightlines = wv_parameters['nbrightlines']
    else:
        nbrightlines = [int(idum) for idum in args.nbrightlines.split(',')]
    poly_crval1_linear = wv_parameters['poly_crval1_linear']
    poly_cdelt1_linear = wv_parameters['poly_cdelt1_linear']

    # list of slitlets to be computed
    list_slitlets = range(islitlet_min, islitlet_max + 1)

    # read master arc line wavelengths (only brightest lines)
    wv_master = read_wv_master_file(
        wv_master_file=args.wv_master_file,
        lines='brightest',
        debugplot=args.debugplot
    )

    # read master arc line wavelengths (whole data set)
    wv_master_all = read_wv_master_file(
        wv_master_file=args.wv_master_file,
        lines='all',
        debugplot=args.debugplot
    )

    # check that the arc lines in the master file are properly sorted
    # in ascending order
    for i in range(len(wv_master_all) - 1):
        if wv_master_all[i] >= wv_master_all[i+1]:
            print('>>> wavelengths: ', wv_master_all[i], wv_master_all[i+1])
            raise ValueError("Arc lines are not properly sorted in master "
                             "file:\n" + args.wv_master_file)

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
            ymargin_bb=args.ymargin_bb,
            debugplot=args.debugplot
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
        if args.remove_sp_background:
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
            times_sigma_threshold=args.times_sigma_threshold)

        # continue working with current slitlet only if arc lines have
        # been detected
        if slt.list_arc_lines is not None:

            # compute intersections between spectrum trails and arc lines
            slt.xy_spectrail_arc_intersections(slitlet2d=slitlet2d)

            # compute rectification transformation
            slt.estimate_tt_to_rectify(order=args.order_fmap,
                                       slitlet2d=slitlet2d)

            # rectify image
            slitlet2d_rect = slt.rectify(slitlet2d,
                                         resampling=2,
                                         transformation=1)

            # median spectrum and line peaks from rectified image
            sp_median, fxpeaks = slt.median_spectrum_from_rectified_image(
                slitlet2d_rect,
                sigma_gaussian_filtering=args.sigma_gaussian_filtering,
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
            expected_wvmin = crval1_linear - args.margin_npix * cdelt1_linear
            naxis1_linear = sp_median.shape[0]
            crvaln_linear = crval1_linear + (naxis1_linear - 1) * cdelt1_linear
            expected_wvmax = crvaln_linear + args.margin_npix * cdelt1_linear

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
                poly_degree_wfit=args.poldeg_initial,
                wv_master=wv_master_eff,
                wv_ini_search=expected_wvmin,
                wv_end_search=expected_wvmax,
                geometry=geometry,
                debugplot=slt.debugplot
            )
            # store initial wavelength calibration polynomial in current
            # slitlet instance
            slt.wpoly = np.polynomial.Polynomial(solution_wv.coeff)
            pause_debugplot(args.debugplot)

            # clip initial master arc line list with all the lines to expected
            # wavelength range
            lok1 = expected_wvmin <= wv_master_all
            lok2 = wv_master_all <= expected_wvmax
            lok = lok1 * lok2
            wv_master_all_eff = wv_master_all[lok]

            # refine wavelength calibration
            if args.poldeg_refined > 0:
                plottitle = args.fitsfile.name + \
                            ' [slitlet#{}, refined]'.format(islitlet)
                poly_refined, yres_summary = refine_arccalibration(
                    sp=sp_median,
                    poly_initial=slt.wpoly,
                    wv_master=wv_master_all_eff,
                    poldeg=args.poldeg_refined,
                    ntimes_match_wv=1,
                    interactive=args.interactive,
                    threshold=args.threshold_wv,
                    plottitle=plottitle,
                    ylogscale=args.ylogscale,
                    geometry=geometry,
                    pdf=pdf,
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

        if args.debugplot == 0:
            if islitlet % 10 == 0:
                if cout != 'x':
                    cout = str(islitlet // 10)
            sys.stdout.write(cout)
            # print(slt)
            if islitlet == list_slitlets[-1]:
                sys.stdout.write('\n')
            sys.stdout.flush()
        else:
            pause_debugplot(args.debugplot)

    # ---

    # Save image with collapsed spectra employed to determine the
    # wavelength calibration
    if args.out_55sp is not None:
        save_ndarray_to_fits(
            array=image2d_55sp,
            file_name=args.out_55sp,
            cast_to_float=True,
            overwrite=True
        )

    # ---

    # Generate structure to save results into an ouptut JSON file
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
    rectwv_coeff.writeto(args.out_json.name)
    print('>>> Saving file ' + args.out_json.name)
    # debugging __getstate__ and __setstate__
    # check_setstate_getstate(rectwv_coeff, args.out_json.name)

    # ---

    if args.out_rect is not None:

        image2d_rectified = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))
        image2d_unrectified = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

        for slt in measured_slitlets:

            islitlet = slt.islitlet

            if abs(args.debugplot) >= 10:
                print(slt)
            else:
                if islitlet % 10 == 0:
                    cout = str(islitlet // 10)
                else:
                    cout = '.'
                sys.stdout.write(cout)
                if islitlet == list_slitlets[-1]:
                    sys.stdout.write('\n')
                sys.stdout.flush()

            # minimum and maximum useful scan (pixel in the spatial direction)
            # for the rectified slitlet
            nscan_min, nscan_max = nscan_minmax_frontiers(
                slt.y0_frontier_lower,
                slt.y0_frontier_upper,
                resize=False
            )

            # extract 2D image corresponding to the selected slitlet: note that
            # in this case we are not using select_unrectified_slitlets()
            # because it introduces extra zero pixels in the slitlet frontiers
            slitlet2d = slt.extract_slitlet2d(image2d)

            # rectify image
            if slt.ttd_order is not None:
                transformation = 1
            elif slt.ttd_order_longslit_model is not None:
                transformation = 2
            else:
                raise ValueError("No ttd transformation defined!")

            slitlet2d_rect = slt.rectify(slitlet2d,
                                         resampling=2,
                                         transformation=transformation)

            slitlet2d_unrect = slt.rectify(slitlet2d_rect,
                                           resampling=2,
                                           transformation=transformation,
                                           inverse=True)

            ii1 = nscan_min - slt.bb_ns1_orig
            ii2 = nscan_max - slt.bb_ns1_orig + 1

            j1 = slt.bb_nc1_orig - 1
            j2 = slt.bb_nc2_orig
            i1 = slt.bb_ns1_orig - 1 + ii1
            i2 = i1 + ii2 - ii1

            image2d_rectified[i1:i2, j1:j2] = slitlet2d_rect[ii1:ii2, :]
            image2d_unrectified[i1:i2, j1:j2] = slitlet2d_unrect[ii1:ii2, :]

        if abs(args.debugplot) % 10 != 0:
            ximshow(image2d_rectified, debugplot=12)
            ximshow(image2d_unrectified, debugplot=12)

        # Save rectified (but not wavelength calibrated image) in first
        # extension, while the second extension is employed to store
        # the unrectified version of the previous one
        save_ndarray_to_fits(
            array=[image2d_rectified, image2d_unrectified],
            file_name=args.out_rect,
            cast_to_float=[True] * 2,
            overwrite=True
        )

    if pdf is not None:
        pdf.close()


if __name__ == "__main__":
    main()
