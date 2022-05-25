#
# Copyright 2018-2020 Universidad Complutense de Madrid
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

"""
Rectification and wavelength calibration polynomials from empirical library
"""

from __future__ import division, print_function

import argparse
from astropy.io import fits
from datetime import datetime
import logging
import numpy as np
from scipy.interpolate import interp1d
import sys
from uuid import uuid4

from numina.array.distortion import ncoef_fmap
from numina.tools.arg_file_is_new import arg_file_is_new
import numina.types.qc

import emirdrp.datamodel as datamodel
from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.instrument.dtuconf import DtuConf
from emirdrp.products import MasterRectWave
from emirdrp.products import RectWaveCoeff
from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.fit_boundaries import expected_distorted_boundaries
from emirdrp.tools.fit_boundaries import expected_distorted_frontiers
from emirdrp.processing.wavecal.slitlet2darc import expected_y0_lower_frontier
from emirdrp.processing.wavecal.slitlet2darc import expected_y0_upper_frontier

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS


def rectwv_coeff_from_mos_library(reduced_image,
                                  master_rectwv,
                                  ignore_dtu_configuration=True,
                                  debugplot=0):
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
    ignore_dtu_configuration : bool
        If True, ignore differences in DTU configuration.
    debugplot : int
        Debugging level for messages and plots. For details see
        'numina.array.display.pause_debugplot.py'.

    Returns
    -------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.

    """

    logger = logging.getLogger(__name__)
    logger.info('Computing expected RectWaveCoeff from CSU configuration')

    # header
    header = reduced_image[0].header
    mecs_header = datamodel.get_mecs_header(reduced_image)

    # read the CSU configuration from the image header
    csu_conf = CsuConfiguration.define_from_header(mecs_header)

    # read the DTU configuration from the image header
    dtu_conf = DtuConf.from_img(reduced_image)

    # retrieve DTU configuration from MasterRectWave object
    dtu_conf_calib = DtuConf.from_values(
        **master_rectwv.meta_info['dtu_configuration']
    )
    # check that the DTU configuration employed to obtain the calibration
    # corresponds to the DTU configuration in the input FITS file
    if dtu_conf != dtu_conf_calib:
        if ignore_dtu_configuration:
            logger.warning('DTU configuration differences found!')
        else:
            logger.info('DTU configuration from image header:')
            logger.info(dtu_conf)
            logger.info('DTU configuration from master calibration:')
            logger.info(dtu_conf_calib)
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
    outdict['meta_info'] = {}
    outdict['meta_info']['creation_date'] = datetime.now().isoformat()
    outdict['meta_info']['description'] = \
        'computation of rectification and wavelength calibration polynomial ' \
        'coefficients for a particular CSU configuration from a MOS model '
    outdict['meta_info']['recipe_name'] = 'undefined'
    outdict['meta_info']['origin'] = {}
    outdict['meta_info']['origin']['fits_frame_uuid'] = 'TBD'
    outdict['meta_info']['origin']['rect_wpoly_mos_uuid'] = \
        master_rectwv.uuid
    outdict['meta_info']['origin']['fitted_boundary_param_uuid'] = \
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
    if abs(debugplot) >= 10:
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
        # employed in Slitlet2dArc.__init__()
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

    # additional parameters (see Slitlet2dArc.__init__)
    for islitlet in list_valid_islitlets:
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        # define expected frontier ordinates at x0_reference for the rectified
        # image imposing the vertical length of the slitlet to be constant
        # and equal to EMIR_NPIXPERSLIT_RECTIFIED
        outdict['contents'][cslitlet]['y0_frontier_lower_expected'] = \
            expected_y0_lower_frontier(islitlet)
        outdict['contents'][cslitlet]['y0_frontier_upper_expected'] = \
            expected_y0_upper_frontier(islitlet)
        # compute linear transformation to place the rectified slitlet at
        # the center of the current slitlet bounding box
        tmpdict = outdict['contents'][cslitlet]
        xdum1 = tmpdict['y0_frontier_lower']
        ydum1 = tmpdict['y0_frontier_lower_expected']
        xdum2 = tmpdict['y0_frontier_upper']
        ydum2 = tmpdict['y0_frontier_upper_expected']
        corr_yrect_b = (ydum2 - ydum1) / (xdum2 - xdum1)
        corr_yrect_a = ydum1 - corr_yrect_b * xdum1
        # compute expected location of rectified boundaries
        y0_reference_lower_expected = \
            corr_yrect_a + corr_yrect_b * tmpdict['y0_reference_lower']
        y0_reference_middle_expected = \
            corr_yrect_a + corr_yrect_b * tmpdict['y0_reference_middle']
        y0_reference_upper_expected = \
            corr_yrect_a + corr_yrect_b * tmpdict['y0_reference_upper']
        # shift transformation to center the rectified slitlet within the
        # slitlet bounding box
        ydummid = (ydum1 + ydum2) / 2
        ioffset = int(
            ydummid - (tmpdict['bb_ns1_orig'] + tmpdict['bb_ns2_orig']) / 2.0)
        corr_yrect_a -= ioffset
        # minimum and maximum row in the rectified slitlet encompassing
        # EMIR_NPIXPERSLIT_RECTIFIED pixels
        # a) scan number (in pixels, from 1 to NAXIS2)
        xdum1 = corr_yrect_a + \
                corr_yrect_b * tmpdict['y0_frontier_lower']
        xdum2 = corr_yrect_a + \
                corr_yrect_b * tmpdict['y0_frontier_upper']
        # b) row number (starting from zero)
        min_row_rectified = \
            int((round(xdum1 * 10) + 5) / 10) - tmpdict['bb_ns1_orig']
        max_row_rectified = \
            int((round(xdum2 * 10) - 5) / 10) - tmpdict['bb_ns1_orig']
        # save previous results in outdict
        outdict['contents'][cslitlet]['y0_reference_lower_expected'] = \
            y0_reference_lower_expected
        outdict['contents'][cslitlet]['y0_reference_middle_expected'] = \
            y0_reference_middle_expected
        outdict['contents'][cslitlet]['y0_reference_upper_expected'] = \
            y0_reference_upper_expected
        outdict['contents'][cslitlet]['corr_yrect_a'] = corr_yrect_a
        outdict['contents'][cslitlet]['corr_yrect_b'] = corr_yrect_b
        outdict['contents'][cslitlet]['min_row_rectified'] = min_row_rectified
        outdict['contents'][cslitlet]['max_row_rectified'] = max_row_rectified

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
                'csu_bar_left': csu_conf.csu_bar_left(islitlet),
                'csu_bar_right': csu_conf.csu_bar_right(islitlet),
                'csu_bar_slit_center': csu_conf.csu_bar_slit_center(islitlet),
                'csu_bar_slit_width': csu_conf.csu_bar_slit_width(islitlet),
                'x0_reference': float(EMIR_NAXIS1) / 2.0 + 0.5,
                'y0_frontier_lower_expected':
                    expected_y0_lower_frontier(islitlet),
                'y0_frontier_upper_expected':
                    expected_y0_upper_frontier(islitlet)
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


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: evaluate rectification and wavelength '
                    'calibration polynomials for the CSU configuration of a '
                    'particular image'
    )

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file",
                        type=argparse.FileType('rb'))
    parser.add_argument("--rect_wpoly_MOSlibrary", required=True,
                        help="Input JSON file with library of rectification "
                             "and wavelength calibration coefficients",
                        type=argparse.FileType('rt'))
    parser.add_argument("--out_json", required=True,
                        help="Output JSON file with calibration computed for "
                             "the input FITS file",
                        type=lambda x: arg_file_is_new(parser, x, mode='wt'))
    # optional arguments
    parser.add_argument("--global_integer_offset_x_pix",
                        help="Global integer offset in the X direction "
                             "(default=0)",
                        default=0, type=int)
    parser.add_argument("--global_integer_offset_y_pix",
                        help="Global integer offset in the Y direction "
                             "(default=0)",
                        default=0, type=int)
    parser.add_argument("--ignore_dtu_configuration",
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

    # ---

    # generate HDUList object
    hdulist = fits.open(args.fitsfile)

    # generate MasterRectWave object
    master_rectwv = MasterRectWave._datatype_load(
        args.rect_wpoly_MOSlibrary.name)

    # compute rectification and wavelength calibration coefficients
    rectwv_coeff = rectwv_coeff_from_mos_library(
        hdulist,
        master_rectwv,
        ignore_dtu_configuration=args.ignore_dtu_configuration,
        debugplot=args.debugplot
    )

    # set global offsets
    rectwv_coeff.global_integer_offset_x_pix = \
        args.global_integer_offset_x_pix
    rectwv_coeff.global_integer_offset_y_pix = \
        args.global_integer_offset_y_pix

    # save RectWaveCoeff object into JSON file
    rectwv_coeff.writeto(args.out_json.name)
    print('>>> Saving file ' + args.out_json.name)
    # debugging __getstate__ and __setstate__
    # check_setstate_getstate(rectwv_coeff, args.out_json.name)


if __name__ == "__main__":
    main()
