from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
from datetime import datetime
import json
import numpy as np
from scipy.interpolate import interp1d
import sys
from uuid import uuid4

from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.instrument.dtu_configuration import DtuConfiguration

from arg_file_is_new import arg_file_is_new
from fit_boundaries import bound_params_from_dict
from fit_boundaries import expected_distorted_boundaries
from fit_boundaries import expected_distorted_frontiers
from rect_wpoly_for_mos import islitlet_progress

from numina.array.distortion import ncoef_fmap
from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from numina.array.distortion import order_fmap
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(prog='evaluate_rect_wpoly')
    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file",
                        type=argparse.FileType('r'))
    parser.add_argument("--rect_wpoly_MOSlibrary", required=True,
                        help="Input JSON file with library of rectification "
                             "and wavelength calibration coefficients",
                        type=argparse.FileType('r'))
    parser.add_argument("--out_rect_wpoly", required=True,
                        help="Output JSON file with calibration computed for "
                             "the input FITS file",
                        type=lambda x: arg_file_is_new(parser, x))
    # optional arguments
    parser.add_argument("--ignore_DTUconf",
                        help="Ignore DTU configurations differences between "
                             "model and input image",
                        action="store_true")
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

    # read the CSU configuration from the header of the input FITS file
    csu_conf = CsuConfiguration.define_from_fits(args.fitsfile)
    if abs(args.debugplot) >= 10:
        print(csu_conf)

    # read the DTU configuration from the header of the input FITS file
    dtu_conf = DtuConfiguration.define_from_fits(args.fitsfile)

    # read calibration structure from JSON file
    rect_wpoly_dict = json.loads(open(args.rect_wpoly_MOSlibrary.name).read())

    dtu_conf_calib = DtuConfiguration.define_from_dictionary(
        rect_wpoly_dict['dtu_configuration'])
    # check that the DTU configuration employed to obtain the calibration
    # corresponds to the DTU configuration in the input FITS file
    if dtu_conf != dtu_conf_calib:
        print('>>> DTU configuration from FITS header:')
        print(dtu_conf)
        print('>>> DTU configuration from calibration JSON file:')
        print(dtu_conf_calib)
        if args.ignore_DTUconf:
            print('WARNING: DTU configuration differences found!')
        else:
            raise ValueError("DTU configurations do not match!")
    else:
        if abs(args.debugplot) >= 10:
            print('>>> DTU Configuration match!')
            print(dtu_conf)

    # read FITS image and its corresponding header
    hdulist = fits.open(args.fitsfile)
    header = hdulist[0].header
    image2d = hdulist[0].data
    hdulist.close()

    # protections
    naxis2, naxis1 = image2d.shape
    if naxis1 != header['naxis1'] or naxis2 != header['naxis2']:
        print('Something is wrong with NAXIS1 and/or NAXIS2')
    if abs(args.debugplot) >= 10:
        print('>>> NAXIS1:', naxis1)
        print('>>> NAXIS2:', naxis2)

    # check that the input FITS file grism and filter match
    filter_name = header['filter']
    if filter_name != rect_wpoly_dict['tags']['filter']:
        raise ValueError("Filter name does not match!")
    grism_name = header['grism']
    if grism_name != rect_wpoly_dict['tags']['grism']:
        raise ValueError("Filter name does not match!")
    if abs(args.debugplot) >= 10:
        print('>>> grism.......:', grism_name)
        print('>>> filter......:', filter_name)

    # read islitlet_min and islitlet_max from input JSON file
    islitlet_min = rect_wpoly_dict['tags']['islitlet_min']
    islitlet_max = rect_wpoly_dict['tags']['islitlet_max']
    if abs(args.debugplot) >= 10:
        print('>>> islitlet_min:', islitlet_min)
        print('>>> islitlet_max:', islitlet_max)

    # Initialize structure to save results into an ouptut JSON file
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
        rect_wpoly_dict['uuid']
    outdict['meta-info']['origin']['fitted_boundary_param_uuid'] = \
        rect_wpoly_dict['fitted_bound_param']['uuid']
    outdict['tags'] = {}
    outdict['tags']['grism'] = grism_name
    outdict['tags']['filter'] = filter_name
    outdict['tags']['islitlet_min'] = islitlet_min
    outdict['tags']['islitlet_max'] = islitlet_max
    outdict['dtu_configuration'] = dtu_conf.outdict()
    outdict['uuid'] = str(uuid4())
    outdict['contents'] = {}

    # compute rectification and wavelength calibration coefficients for each
    # slitlet according to its csu_bar_slit_center value
    for islitlet in range(islitlet_min, islitlet_max + 1):
        islitlet_progress(islitlet, islitlet_max)
        cslitlet = 'slitlet' + str(islitlet).zfill(2)

        # csu_bar_slit_center of current slitlet in initial FITS image
        csu_bar_slit_center = csu_conf.csu_bar_slit_center(islitlet)

        # input data structure
        tmpdict = rect_wpoly_dict['contents'][cslitlet]
        list_csu_bar_slit_center = tmpdict['list_csu_bar_slit_center']

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
                                          kind='linear')
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
                                      kind='linear')
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
    fitted_bound_param = rect_wpoly_dict['fitted_bound_param']
    parmodel = fitted_bound_param['meta-info']['parmodel']
    params = bound_params_from_dict(fitted_bound_param)
    if abs(args.debugplot) >= 10:
        print('* Fitted boundary parameters:')
        params.pretty_print()
    for islitlet in range(islitlet_min, islitlet_max + 1):
        islitlet_progress(islitlet, islitlet_max)
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
    for islitlet in range(islitlet_min, islitlet_max + 1):
        islitlet_progress(islitlet, islitlet_max)
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        # parameters already available in the input JSON file
        for par in ['bb_nc1_orig', 'bb_nc2_orig', 'ymargin_bb']:
            outdict['contents'][cslitlet][par] = \
                rect_wpoly_dict['contents'][cslitlet][par]
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
        ymargin_bb = rect_wpoly_dict['contents'][cslitlet]['ymargin_bb']
        bb_ns1_orig = int(ylower.min() + 0.5) - ymargin_bb
        if bb_ns1_orig < 1:
            bb_ns1_orig = 1
        bb_ns2_orig = int(yupper.max() + 0.5) + ymargin_bb
        if bb_ns2_orig > EMIR_NAXIS2:
            bb_ns2_orig = EMIR_NAXIS2
        outdict['contents'][cslitlet]['bb_ns1_orig'] = bb_ns1_orig
        outdict['contents'][cslitlet]['bb_ns2_orig'] = bb_ns2_orig

    # Save resulting JSON structure
    with open(args.out_rect_wpoly.name, 'w') as fstream:
        json.dump(outdict, fstream, indent=2, sort_keys=True)
        print('>>> Saving file ' + args.out_rect_wpoly.name)


if __name__ == "__main__":
    main()
