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
import json
import sys

from numina.tools.arg_file_is_new import arg_file_is_new
from numina.tools.check_setstate_getstate import check_setstate_getstate
import numina.types.qc
from emirdrp.products import RefinedBoundaryModelParam


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: convert JSON file with refined multislit '
                    'parameters to new JSON format'
    )

    # required arguments
    parser.add_argument("input_json",
                        help="Input JSON with refined boundary parameters",
                        type=argparse.FileType('rt'))
    parser.add_argument("output_json",
                        help="Output JSON with fitted boundary parameters",
                        type=lambda x: arg_file_is_new(parser, x, mode='wt'))

    # optional arguments
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")
    args = parser.parse_args(args)

    if args.echo:
        print('\033[1m\033[31m% ' + ' '.join(sys.argv) + '\033[0m\n')

    # ---

    # read input JSON file
    input_json = json.loads(open(args.input_json.name).read())

    # generate object of type RefinedBoundaryModelParam from input JSON file
    refined_boundary_model = RefinedBoundaryModelParam(instrument='EMIR')
    refined_boundary_model.tags = {
        'grism': input_json['tags']['grism'],
        'filter': input_json['tags']['filter']
    }
    refined_boundary_model.contents = input_json['contents']
    refined_boundary_model.meta_info['dtu_configuration'] = \
        input_json['dtu_configuration']
    refined_boundary_model.meta_info['dtu_configuration_maxdiff'] = \
        input_json['dtu_configuration_maxdiff']
    for item in ['function_evaluations', 'global_residual', 'maxDTUoffset',
                 'numresolution', 'parmodel', 'tolerance']:
        refined_boundary_model.meta_info[item] = input_json['meta_info'][item]
    refined_boundary_model.meta_info['origin']['bounddict'] = \
        'uuid:' + input_json['meta_info']['origin']['bounddict_uuid']
    refined_boundary_model.meta_info['origin']['init_bound_param'] = \
        'uuid:' + input_json['meta_info']['origin']['init_bound_param_uuid']
    refined_boundary_model.quality_control = numina.types.qc.QC.GOOD
    refined_boundary_model.writeto(args.output_json.name)
    # debugging __getstate__ and __setstate__
    # check_setstate_getstate(refined_boundary_model, args.output_json.name)


if __name__ == "__main__":
    main()
