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

from __future__ import division, print_function

import argparse
from datetime import datetime
import logging
import numpy as np
import sys
from uuid import uuid4

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.logging_from_debugplot import logging_from_debugplot
from numina.array.display.polfit_residuals import \
    polfit_residuals_with_sigma_rejection
from numina.array.distortion import ncoef_fmap
from numina.array.stats import summary
from numina.tools.arg_file_is_new import arg_file_is_new
from numina.tools.check_setstate_getstate import check_setstate_getstate

from emirdrp.products import RectWaveCoeff

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NBARS


def rectwv_coeff_add_longslit_model(rectwv_coeff, geometry, debugplot=0):
    """Compute longslit_model coefficients for RectWaveCoeff object.

    Parameters
    ----------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for a
        particular CSU configuration corresponding to a longslit
        observation.
    geometry : TBD
    debugplot : int
        Debugging level for messages and plots. For details see
        'numina.array.display.pause_debugplot.py'.

    Returns
    -------
    rectwv_coeff : RectWaveCoeff instance
        Updated object with longslit_model coefficients computed.

    """

    logger = logging.getLogger(__name__)

    # check grism and filter
    grism_name = rectwv_coeff.tags['grism']
    logger.info('Grism: ' + grism_name)
    filter_name = rectwv_coeff.tags['filter']
    logger.info('Filter: ' + filter_name)

    # list of slitlets to be computed
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in rectwv_coeff.missing_slitlets:
        list_valid_islitlets.remove(idel)
    if abs(debugplot) >= 10:
        print('>>> valid slitlet numbers:\n', list_valid_islitlets)

    # ---

    # check that the CSU configuration corresponds to longslit
    csu_bar_slit_center_list = []
    for islitlet in list_valid_islitlets:
        csu_bar_slit_center_list.append(
            rectwv_coeff.contents[islitlet - 1]['csu_bar_slit_center']
        )
    if abs(debugplot) >= 10:
        logger.debug('Checking csu_bar_slit_center values:')
        summary(np.array(csu_bar_slit_center_list), debug=True)
        pause_debugplot(debugplot)

    # ---

    # polynomial coefficients corresponding to the wavelength calibration

    # step 0: determine poldeg_refined, checking that it is the same for
    # all the slitlets
    poldeg_refined_list = []
    for islitlet in list_valid_islitlets:
        poldeg_refined_list.append(
            len(rectwv_coeff.contents[islitlet - 1]['wpoly_coeff']) - 1
        )
    # remove duplicates
    poldeg_refined_list = list(set(poldeg_refined_list))
    if len(poldeg_refined_list) != 1:
        raise ValueError('Unexpected different poldeg_refined found: ' +
                         str(poldeg_refined_list))
    poldeg_refined = poldeg_refined_list[0]

    # step 1: compute variation of each coefficient as a function of
    # y0_reference_middle of each slitlet
    list_poly = []
    for i in range(poldeg_refined + 1):
        xp = []
        yp = []
        for islitlet in list_valid_islitlets:
            tmp_dict = rectwv_coeff.contents[islitlet - 1]
            wpoly_coeff = tmp_dict['wpoly_coeff']
            if wpoly_coeff is not None:
                xp.append(tmp_dict['y0_reference_middle'])
                yp.append(wpoly_coeff[i])
        poly, yres, reject = polfit_residuals_with_sigma_rejection(
            x=np.array(xp),
            y=np.array(yp),
            deg=2,
            times_sigma_reject=5,
            xlabel='y0_rectified',
            ylabel='coeff[' + str(i) + ']',
            title="Fit to refined wavelength calibration coefficients",
            geometry=geometry,
            debugplot=debugplot
        )
        list_poly.append(poly)

    # step 2: use the variation of each polynomial coefficient with
    # y0_reference_middle to infer the expected wavelength calibration
    # polynomial for each rectifified slitlet
    for islitlet in list_valid_islitlets:
        tmp_dict = rectwv_coeff.contents[islitlet - 1]
        y0_reference_middle = tmp_dict['y0_reference_middle']
        list_new_coeff = []
        for i in range(poldeg_refined + 1):
            new_coeff = list_poly[i](y0_reference_middle)
            list_new_coeff.append(new_coeff)
        tmp_dict['wpoly_coeff_longslit_model'] = list_new_coeff

    # ---

    # rectification transformation coefficients aij and bij

    # step 0: determine order_fmap, checking that it is the same for
    # all the slitlets
    order_fmap_list = []
    for islitlet in list_valid_islitlets:
        order_fmap_list.append(
            rectwv_coeff.contents[islitlet - 1]['ttd_order']
        )
    # remove duplicates
    order_fmap_list = list(set(order_fmap_list))
    if len(order_fmap_list) != 1:
        raise ValueError('Unexpected different order_fmap found')
    order_fmap = order_fmap_list[0]

    # step 1: compute variation of each coefficient as a function of
    # y0_reference_middle of each slitlet
    list_poly_ttd_aij = []
    list_poly_ttd_bij = []
    list_poly_tti_aij = []
    list_poly_tti_bij = []
    ncoef_ttd = ncoef_fmap(order_fmap)
    for i in range(ncoef_ttd):
        xp = []
        yp_ttd_aij = []
        yp_ttd_bij = []
        yp_tti_aij = []
        yp_tti_bij = []
        for islitlet in list_valid_islitlets:
            tmp_dict = rectwv_coeff.contents[islitlet - 1]
            ttd_aij = tmp_dict['ttd_aij']
            ttd_bij = tmp_dict['ttd_bij']
            tti_aij = tmp_dict['tti_aij']
            tti_bij = tmp_dict['tti_bij']
            if ttd_aij is not None:
                xp.append(tmp_dict['y0_reference_middle'])
                yp_ttd_aij.append(ttd_aij[i])
                yp_ttd_bij.append(ttd_bij[i])
                yp_tti_aij.append(tti_aij[i])
                yp_tti_bij.append(tti_bij[i])
        poly, yres, reject = polfit_residuals_with_sigma_rejection(
            x=np.array(xp),
            y=np.array(yp_ttd_aij),
            deg=5,
            times_sigma_reject=5,
            xlabel='y0_rectified',
            ylabel='ttd_aij[' + str(i) + ']',
            geometry=geometry,
            debugplot=debugplot
        )
        list_poly_ttd_aij.append(poly)
        poly, yres, reject = polfit_residuals_with_sigma_rejection(
            x=np.array(xp),
            y=np.array(yp_ttd_bij),
            deg=5,
            times_sigma_reject=5,
            xlabel='y0_rectified',
            ylabel='ttd_bij[' + str(i) + ']',
            geometry=geometry,
            debugplot=debugplot
        )
        list_poly_ttd_bij.append(poly)
        poly, yres, reject = polfit_residuals_with_sigma_rejection(
            x=np.array(xp),
            y=np.array(yp_tti_aij),
            deg=5,
            times_sigma_reject=5,
            xlabel='y0_rectified',
            ylabel='tti_aij[' + str(i) + ']',
            geometry=geometry,
            debugplot=debugplot
        )
        list_poly_tti_aij.append(poly)
        poly, yres, reject = polfit_residuals_with_sigma_rejection(
            x=np.array(xp),
            y=np.array(yp_tti_bij),
            deg=5,
            times_sigma_reject=5,
            xlabel='y0_rectified',
            ylabel='tti_bij[' + str(i) + ']',
            geometry=geometry,
            debugplot=debugplot
        )
        list_poly_tti_bij.append(poly)

    # step 2: use the variation of each coefficient with y0_reference_middle
    # to infer the expected rectification transformation for each slitlet
    for islitlet in list_valid_islitlets:
        tmp_dict = rectwv_coeff.contents[islitlet - 1]
        y0_reference_middle = tmp_dict['y0_reference_middle']
        tmp_dict['ttd_order_longslit_model'] = order_fmap
        ttd_aij_longslit_model = []
        ttd_bij_longslit_model = []
        tti_aij_longslit_model = []
        tti_bij_longslit_model = []
        for i in range(ncoef_ttd):
            new_coeff = list_poly_ttd_aij[i](y0_reference_middle)
            ttd_aij_longslit_model.append(new_coeff)
            new_coeff = list_poly_ttd_bij[i](y0_reference_middle)
            ttd_bij_longslit_model.append(new_coeff)
            new_coeff = list_poly_tti_aij[i](y0_reference_middle)
            tti_aij_longslit_model.append(new_coeff)
            new_coeff = list_poly_tti_bij[i](y0_reference_middle)
            tti_bij_longslit_model.append(new_coeff)
        tmp_dict['ttd_aij_longslit_model'] = ttd_aij_longslit_model
        tmp_dict['ttd_bij_longslit_model'] = ttd_bij_longslit_model
        tmp_dict['tti_aij_longslit_model'] = tti_aij_longslit_model
        tmp_dict['tti_bij_longslit_model'] = tti_bij_longslit_model

    # ---

    # update uuid and meta_info in output JSON structure
    rectwv_coeff.uuid = str(uuid4())
    rectwv_coeff.meta_info['creation_date'] = datetime.now().isoformat()

    # return updated object
    return rectwv_coeff


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument("--input_rectwv_coeff", required=True,
                        help="Input JSON file with rectification and "
                             "wavelength calibration polynomials "
                             "corresponding to a longslit observation",
                        type=argparse.FileType('rt'))
    parser.add_argument("--output_rectwv_coeff", required=True,
                        help="Output JSON file with updated longslit_model "
                             "coefficients",
                        type=lambda x: arg_file_is_new(parser, x, mode='wt'))

    # optional arguments
    parser.add_argument("--geometry",
                        help="tuple x,y,dx,dy (default 0,0,640,480)",
                        default="0,0,640,480")
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

    logging_from_debugplot(args.debugplot)
    logger = logging.getLogger(__name__)

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

    # generate RectWaveCoeff object
    rectwv_coeff = RectWaveCoeff._datatype_load(
        args.input_rectwv_coeff.name)

    # update longslit_model parameters
    rectwv_coeff_updated = rectwv_coeff_add_longslit_model(
        rectwv_coeff=rectwv_coeff,
        geometry=geometry,
        debugplot=args.debugplot
    )

    # save updated RectWaveCoeff object into JSON file
    rectwv_coeff_updated.writeto(args.output_rectwv_coeff.name)
    logger.info('>>> Saving file ' + args.output_rectwv_coeff.name)
    # debugging __getstate__ and __setstate__
    # check_setstate_getstate(rectwv_coeff_updated, args.output_rectwv_coeff.name)


if __name__ == "__main__":
    main()
