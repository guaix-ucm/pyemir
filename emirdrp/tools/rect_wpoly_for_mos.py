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
from datetime import datetime
import numpy as np
import sys
from uuid import uuid4

from numina.array.display.fileinfo import list_fileinfo_from_txt
from numina.array.distortion import ncoef_fmap
from numina.tools.arg_file_is_new import arg_file_is_new
from numina.tools.check_setstate_getstate import check_setstate_getstate
import numina.types.qc

from emirdrp.instrument.dtuconf import DtuConf
from emirdrp.products import RefinedBoundaryModelParam
from emirdrp.products import RectWaveCoeff
from emirdrp.products import MasterRectWave

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NBARS


def islitlet_progress(islitlet, islitlet_max, ignore=False):
    """Auxiliary function to print out progress in loop of slitlets.

    Parameters
    ----------
    islitlet : int
        Current slitlet number.
    islitlet_max : int
        Maximum slitlet number.
    ignore : bool
        If True, display 'i'

    """

    if ignore:
        cout = 'i'
    else:
        if islitlet % 10 == 0:
            cout = str(islitlet // 10)
        else:
            cout = '.'

    sys.stdout.write(cout)
    if islitlet == islitlet_max:
        sys.stdout.write('\n')
    sys.stdout.flush()


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(prog='rect_wpoly_for_mos')
    # required arguments
    parser.add_argument("input_list",
                        help="TXT file with list JSON files derived from "
                             "longslit data")
    parser.add_argument("--fitted_bound_param", required=True,
                        help="Input JSON with fitted boundary parameters",
                        type=argparse.FileType('rt'))
    parser.add_argument("--out_MOSlibrary", required=True,
                        help="Output JSON file with results",
                        type=lambda x: arg_file_is_new(parser, x))
    # optional arguments
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

    # Read input TXT file with list of JSON files
    list_json_files = list_fileinfo_from_txt(args.input_list)
    nfiles = len(list_json_files)
    if abs(args.debugplot) >= 10:
        print('>>> Number of input JSON files:', nfiles)
        for item in list_json_files:
            print(item)
    if nfiles < 2:
        raise ValueError("Insufficient number of input JSON files")

    # read fitted boundary parameters and check that all the longslit JSON
    # files have been computed using the same fitted boundary parameters
    refined_boundary_model = RefinedBoundaryModelParam._datatype_load(
        args.fitted_bound_param.name)
    for ifile in range(nfiles):
        coef_rect_wpoly = RectWaveCoeff._datatype_load(
            list_json_files[ifile].filename)
        uuid_tmp = coef_rect_wpoly.meta_info['origin']['bound_param']
        if uuid_tmp[4:] != refined_boundary_model.uuid:
            print('Expected uuid:', refined_boundary_model.uuid)
            print('uuid for ifile #' + str(ifile + 1) + ": " + uuid_tmp)
            raise ValueError("Fitted boundary parameter uuid's do not match")

    # check consistency of grism, filter, DTU configuration and list of
    # valid slitlets
    coef_rect_wpoly_first_longslit = RectWaveCoeff._datatype_load(
        list_json_files[0].filename)
    filter_name = coef_rect_wpoly_first_longslit.tags['filter']
    grism_name = coef_rect_wpoly_first_longslit.tags['grism']
    dtu_conf = DtuConf.from_values(
        **coef_rect_wpoly_first_longslit.meta_info['dtu_configuration']
    )
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in coef_rect_wpoly_first_longslit.missing_slitlets:
        list_valid_islitlets.remove(idel)
    for ifile in range(1, nfiles):
        coef_rect_wpoly = RectWaveCoeff._datatype_load(
            list_json_files[ifile].filename)
        filter_tmp = coef_rect_wpoly.tags['filter']
        if filter_name != filter_tmp:
            print(filter_name)
            print(filter_tmp)
            raise ValueError("Unexpected different filter found")
        grism_tmp = coef_rect_wpoly.tags['grism']
        if grism_name != grism_tmp:
            print(grism_name)
            print(grism_tmp)
            raise ValueError("Unexpected different grism found")
        coef_rect_wpoly = RectWaveCoeff._datatype_load(
            list_json_files[ifile].filename)
        dtu_conf_tmp = DtuConf.from_values(
            **coef_rect_wpoly.meta_info['dtu_configuration']
        )
        if dtu_conf != dtu_conf_tmp:
            print(dtu_conf)
            print(dtu_conf_tmp)
            raise ValueError("Unexpected different DTU configurations found")
        list_valid_islitlets_tmp = list(range(1, EMIR_NBARS + 1))
        for idel in coef_rect_wpoly.missing_slitlets:
            list_valid_islitlets_tmp.remove(idel)
        if list_valid_islitlets != list_valid_islitlets_tmp:
            print(list_valid_islitlets)
            print(list_valid_islitlets_tmp)
            raise ValueError("Unexpected different list of valid slitlets")

    # check consistency of horizontal bounding box limits (bb_nc1_orig and
    # bb_nc2_orig) and ymargin_bb, and store the values for each slitlet
    dict_bb_param = {}
    print("Checking horizontal bounding box limits and ymargin_bb:")
    for islitlet in list(range(1, EMIR_NBARS + 1)):
        if islitlet in list_valid_islitlets:
            islitlet_progress(islitlet, EMIR_NBARS, ignore=False)
            cslitlet = 'slitlet' + str(islitlet).zfill(2)
            dict_bb_param[cslitlet] = {}
            for par in ['bb_nc1_orig', 'bb_nc2_orig', 'ymargin_bb']:
                value_initial = \
                    coef_rect_wpoly_first_longslit.contents[islitlet - 1][par]
                for ifile in range(1, nfiles):
                    coef_rect_wpoly = RectWaveCoeff._datatype_load(
                        list_json_files[ifile].filename)
                    value_tmp = coef_rect_wpoly.contents[islitlet - 1][par]
                    if value_initial != value_tmp:
                        print(islitlet, value_initial, value_tmp)
                        print(value_tmp)
                        raise ValueError("Unexpected different " + par)
                    dict_bb_param[cslitlet][par] = value_initial
        else:
            islitlet_progress(islitlet, EMIR_NBARS, ignore=True)
    print('OK!')

    # ---

    # Read and store all the longslit data
    list_coef_rect_wpoly = []
    for ifile in range(nfiles):
        coef_rect_wpoly = RectWaveCoeff._datatype_load(
            list_json_files[ifile].filename)
        list_coef_rect_wpoly.append(coef_rect_wpoly)

    # ---

    # Initialize structure to save results into an ouptut JSON file
    outdict = {}
    outdict['refined_boundary_model'] = refined_boundary_model.__getstate__()
    outdict['instrument'] = 'EMIR'
    outdict['meta_info'] = {}
    outdict['meta_info']['creation_date'] = datetime.now().isoformat()
    outdict['meta_info']['description'] = \
        'rectification and wavelength calibration polynomial coefficients ' \
        'as a function of csu_bar_slit_center for MOS'
    outdict['meta_info']['recipe_name'] = 'undefined'
    outdict['meta_info']['origin'] = {}
    outdict['meta_info']['origin']['wpoly_longslits'] = {}
    for ifile in range(nfiles):
        cdum = 'longslit_' + str(ifile + 1).zfill(3) + '_uuid'
        outdict['meta_info']['origin']['wpoly_longslits'][cdum] = \
            list_coef_rect_wpoly[ifile].uuid
    outdict['tags'] = {}
    outdict['tags']['grism'] = grism_name
    outdict['tags']['filter'] = filter_name
    outdict['dtu_configuration'] = dtu_conf.outdict()
    outdict['uuid'] = str(uuid4())
    outdict['contents'] = {}

    # include bb_nc1_orig, bb_nc2_orig and ymargin_bb for each slitlet
    # (note that the values of bb_ns1_orig and bb_ns2_orig cannot be
    # computed at this stage because they depend on csu_bar_slit_center)
    for islitlet in list_valid_islitlets:
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        outdict['contents'][cslitlet] = dict_bb_param[cslitlet]

    # check that order for rectification transformations is the same for all
    # the slitlets and longslit configurations
    order_check_list = []
    for ifile in range(nfiles):
        tmpdict = list_coef_rect_wpoly[ifile].contents
        for islitlet in list_valid_islitlets:
            ttd_order = tmpdict[islitlet - 1]['ttd_order']
            if ttd_order is not None:
                order_check_list.append(ttd_order)
            ttd_order_modeled = \
                tmpdict[islitlet - 1]['ttd_order_longslit_model']
            order_check_list.append(ttd_order_modeled)
    # remove duplicates in list
    order_no_duplicates = list(set(order_check_list))
    if len(order_no_duplicates) != 1:
        print('order_no_duplicates:', order_no_duplicates)
        raise ValueError('tdd_order is not constant!')
    ttd_order = int(order_no_duplicates[0])
    ncoef_rect = ncoef_fmap(ttd_order)
    if abs(args.debugplot) >= 10:
        print('>>> ttd_order........:', ttd_order)
        print('>>> ncoef_rect.......:', ncoef_rect)

    # check that polynomial degree in frontiers and spectrails are the same
    poldeg_check_list = []
    for ifile in range(nfiles):
        tmpdict = list_coef_rect_wpoly[ifile].contents
        for islitlet in list_valid_islitlets:
            tmppoly = tmpdict[islitlet - 1]['frontier']['poly_coef_lower']
            poldeg_check_list.append(len(tmppoly) - 1)
            tmppoly = tmpdict[islitlet - 1]['frontier']['poly_coef_upper']
            poldeg_check_list.append(len(tmppoly) - 1)
            tmppoly = tmpdict[islitlet - 1]['spectrail']['poly_coef_lower']
            poldeg_check_list.append(len(tmppoly) - 1)
            tmppoly = tmpdict[islitlet - 1]['spectrail']['poly_coef_middle']
            poldeg_check_list.append(len(tmppoly) - 1)
            tmppoly = tmpdict[islitlet - 1]['spectrail']['poly_coef_upper']
            poldeg_check_list.append(len(tmppoly) - 1)
    # remove duplicates in list
    poldeg_no_duplicates = list(set(poldeg_check_list))
    if len(poldeg_no_duplicates) != 1:
        print('poldeg_no_duplicates:', poldeg_no_duplicates)
        raise ValueError('poldeg is not constant in frontiers and '
                         'spectrails!')
    poldeg_spectrails = int(poldeg_no_duplicates[0])
    if abs(args.debugplot) >= 10:
        print('>>> poldeg spectrails:', poldeg_spectrails)

    # check that polynomial degree of wavelength calibration is the same for
    # all the slitlets
    poldeg_check_list = []
    for ifile in range(nfiles):
        tmpdict = list_coef_rect_wpoly[ifile].contents
        for islitlet in list_valid_islitlets:
            tmppoly = tmpdict[islitlet - 1]['wpoly_coeff']
            poldeg_check_list.append(len(tmppoly) - 1)
            tmppoly = tmpdict[islitlet - 1]['wpoly_coeff_longslit_model']
            poldeg_check_list.append(len(tmppoly) - 1)
    # remove duplicates in list
    poldeg_no_duplicates = list(set(poldeg_check_list))
    if len(poldeg_no_duplicates) != 1:
        print('poldeg_no_duplicates:', poldeg_no_duplicates)
        raise ValueError('poldeg is not constant in wavelength calibration '
                         'polynomials!')
    poldeg_wavecal = int(poldeg_no_duplicates[0])
    if abs(args.debugplot) >= 10:
        print('>>> poldeg wavecal...:', poldeg_wavecal)

    # ---

    # csu_bar_slit_center values for each slitlet
    print("CSU_bar_slit_center values:")
    for islitlet in list(range(1, EMIR_NBARS + 1)):
        if islitlet in list_valid_islitlets:
            islitlet_progress(islitlet, EMIR_NBARS, ignore=False)
            cslitlet = 'slitlet' + str(islitlet).zfill(2)
            list_csu_bar_slit_center = []
            for ifile in range(nfiles):
                tmpdict = list_coef_rect_wpoly[ifile].contents[islitlet - 1]
                csu_bar_slit_center = tmpdict['csu_bar_slit_center']
                list_csu_bar_slit_center.append(csu_bar_slit_center)
            # check that list_csu_bar_slit_center is properly sorted
            if not np.all(list_csu_bar_slit_center[:-1] <=
                      list_csu_bar_slit_center[1:]):
                print('cslitlet: ', cslitlet)
                print('list_csu_bar_slit_center: ', list_csu_bar_slit_center)
                raise ValueError('Unsorted list_csu_bar_slit_center')
            outdict['contents'][cslitlet]['list_csu_bar_slit_center'] = \
                list_csu_bar_slit_center
        else:
            islitlet_progress(islitlet, EMIR_NBARS, ignore=True)
    print('OK!')

    # ---

    # rectification polynomial coefficients

    # note: when aij and bij have not been computed, we use the modeled
    # version aij_longslit_model and bij_longslit_model
    print("Rectification polynomial coefficients:")
    for islitlet in list(range(1, EMIR_NBARS + 1)):
        if islitlet in list_valid_islitlets:
            islitlet_progress(islitlet, EMIR_NBARS, ignore=False)
            cslitlet = 'slitlet' + str(islitlet).zfill(2)
            outdict['contents'][cslitlet]['ttd_order'] = ttd_order
            outdict['contents'][cslitlet]['ncoef_rect'] = ncoef_rect
            for keycoef in ['ttd_aij', 'ttd_bij', 'tti_aij', 'tti_bij']:
                for icoef in range(ncoef_rect):
                    ccoef = str(icoef).zfill(2)
                    list_cij = []
                    for ifile in range(nfiles):
                        tmpdict = \
                            list_coef_rect_wpoly[ifile].contents[islitlet - 1]
                        cij = tmpdict[keycoef]
                        if cij is not None:
                            list_cij.append(cij[icoef])
                        else:
                            cij_modeled = tmpdict[keycoef + '_longslit_model']
                            if cij_modeled is None:
                                raise ValueError("Unexpected cij_modeled=None!")
                            else:
                                list_cij.append(cij_modeled[icoef])
                            if abs(args.debugplot) >= 10:
                                print("Warning: using " + keycoef +
                                      "_longslit_model for " + cslitlet +
                                      " in file " +
                                      list_json_files[ifile].filename)
                    cdum = 'list_' + keycoef + '_' + ccoef
                    outdict['contents'][cslitlet][cdum] = list_cij
        else:
            islitlet_progress(islitlet, EMIR_NBARS, ignore=True)

    print('OK!')

    # ---

    # wavelength calibration polynomial coefficients

    # note: when wpoly_coeff have not been computed, we use the
    # wpoly_coeff_longslit_model
    print("Wavelength calibration polynomial coefficients:")
    for islitlet in list(range(1, EMIR_NBARS + 1)):
        if islitlet in list_valid_islitlets:
            islitlet_progress(islitlet, EMIR_NBARS, ignore=False)
            cslitlet = 'slitlet' + str(islitlet).zfill(2)
            outdict['contents'][cslitlet]['wpoly_degree'] = poldeg_wavecal
            for icoef in range(poldeg_wavecal + 1):
                ccoef = str(icoef).zfill(2)
                list_cij = []
                for ifile in range(nfiles):
                    tmpdict = list_coef_rect_wpoly[ifile].contents[islitlet - 1]
                    cij = tmpdict['wpoly_coeff']
                    if cij is not None:
                        list_cij.append(cij[icoef])
                    else:
                        cij_modeled = tmpdict['wpoly_coeff_longslit_model']
                        if cij_modeled is None:
                            raise ValueError("Unexpected cij_modeled=None!")
                        else:
                            list_cij.append(cij_modeled[icoef])
                        if abs(args.debugplot) >= 10:
                            print("Warning: using wpoly_coeff_longslit_model" +
                                  " for " + cslitlet +
                                  " in file " +
                                  list_json_files[ifile].filename)
                outdict['contents'][cslitlet]['list_wpoly_coeff_' + ccoef] = \
                    list_cij
        else:
            islitlet_progress(islitlet, EMIR_NBARS, ignore=True)
    print('OK!')

    # ---

    # OBSOLETE
    # Save resulting JSON structure
    '''
    with open(args.out_MOSlibrary.name + '_old', 'w') as fstream:
        json.dump(outdict, fstream, indent=2, sort_keys=True)
        print('>>> Saving file ' + args.out_MOSlibrary.name + '_old')
    '''

    # --

    # Create object of type MasterRectWave with library of coefficients
    # for rectification and wavelength calibration
    master_rectwv = MasterRectWave(instrument='EMIR')
    master_rectwv.quality_control = numina.types.qc.QC.GOOD
    master_rectwv.tags['grism'] = grism_name
    master_rectwv.tags['filter'] = filter_name
    master_rectwv.meta_info['dtu_configuration'] = outdict['dtu_configuration']
    master_rectwv.meta_info['refined_boundary_model'] = {
        'parmodel': refined_boundary_model.meta_info['parmodel']
    }
    master_rectwv.meta_info['refined_boundary_model'].update(
        outdict['refined_boundary_model']['contents']
    )
    master_rectwv.total_slitlets = EMIR_NBARS
    master_rectwv.meta_info['origin'] = {
        'bound_param': 'uuid' + refined_boundary_model.uuid,
        'longslit_frames': ['uuid:' + list_coef_rect_wpoly[ifile].uuid
                            for ifile in range(nfiles)]
    }
    for i in range(EMIR_NBARS):
        islitlet = i + 1
        dumdict = {'islitlet': islitlet}
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        if cslitlet in outdict['contents']:
            dumdict.update(outdict['contents'][cslitlet])
        else:
            dumdict.update({
                'bb_nc1_orig': 0,
                'bb_nc2_orig': 0,
                'ymargin_bb': 0,
                'list_csu_bar_slit_center': [],
                'ttd_order': 0,
                'ncoef_rect': 0,
                'wpolydegree': 0
            })
            master_rectwv.missing_slitlets.append(islitlet)
        master_rectwv.contents.append(dumdict)
    master_rectwv.writeto(args.out_MOSlibrary.name)
    print('>>> Saving file ' + args.out_MOSlibrary.name)
    # debugging __getstate__ and __setstate__
    # check_setstate_getstate(master_rectwv, args.out_MOSlibrary.name)


if __name__ == "__main__":
    main()
