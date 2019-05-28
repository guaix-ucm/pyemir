#
# Copyright 2019 Universidad Complutense de Madrid
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

"""Generate YAML files to reduce ABBA observations"""

import argparse
import sys
import textwrap

from numina.array.display.fileinfo import FileInfo
from numina.tools.arg_file_is_new import arg_file_is_new


# auxiliary function to generate content of YAML files
def generate_yaml_content(step_number, args, list_fileinfo, enabled=True):

    # obsid_prefix
    if args.obsid_prefix is None:
        obsid_prefix = ''
    else:
        obsid_prefix = args.obsid_prefix + '_'

    nimages = len(list_fileinfo)

    output = ''
    if step_number == 0:
        # initial rectification and wavelength calibration coefficients
        idlabel = list_fileinfo[0].filename[:10]
        if obsid_prefix == '':
            output += 'id: _{}_rectwv_preliminary\n'.format(idlabel)
        else:
            output += 'id: _{}_{}_rectwv_preliminary\n'.format(obsid_prefix,
                                                               idlabel)
        output += 'instrument: EMIR\n'
        output += 'mode: GENERATE_RECTWV_COEFF\n'
        output += 'frames:\n'
        for i in range(args.npreliminary):
            output += ' - ' + list_fileinfo[i].filename + '\n'
        if enabled:
            output += 'enabled: True'
        else:
            output += 'enabled: False'
    elif step_number == 1:
        # check that all the expected parameters have been provided
        lrequirements = [
            'refine_wavecalib_mode',
            'minimum_slitlet_width_mm',
            'maximum_slitlet_width_mm',
            'global_integer_offset_x_pix',
            'global_integer_offset_y_pix'
        ]
        for item in lrequirements:
            if args.__dict__[item] is None:
                raise ValueError('Parameter {} is None!'.format(item))
        # refined rectification and wavelength calibration for each block
        for i in range(nimages):
            idlabel = list_fileinfo[i].filename[:10]
            if obsid_prefix == '':
                output += 'id: _{}_rectwv\n'.format(idlabel)
            else:
                output += 'id: _{}_{}_rectwv\n'.format(obsid_prefix, idlabel)
            output += 'instrument: EMIR\n'
            output += 'mode: GENERATE_RECTWV_COEFF\n'
            output += 'frames:\n'
            output += ' - ' + list_fileinfo[i].filename + '\n'
            output += 'requirements:\n'
            for item in lrequirements:
                output += '  {}: {}\n'.format(item, args.__dict__[item])
            if enabled:
                output += 'enabled: True'
            else:
                output += 'enabled: False'
            if i < nimages - 1:
                output += '\n---\n'
    elif step_number == 2:
        # apply rectification and wavelength calibration to each block
        output += 'id: _' + obsid_prefix + 'abba_combined\n'
        output += 'instrument: EMIR\n'
        output += 'mode: ABBA_SPECTRA_RECTWV\n'
        output += 'frames:\n'
        for i in range(nimages):
            idlabel = list_fileinfo[i].filename[:10]
            if obsid_prefix == '':
                dumid = '_{}_rectwv'.format(idlabel)
            else:
                dumid = '_{}_{}_rectwv'.format(obsid_prefix, idlabel)
            output += ' - ' + dumid + '\n'
        output += 'requirements:\n'
        output += '  method: sigmaclip\n'
        output += '  refine_objects_in_slit: 1\n'
        output += 'enabled: True\n'
    else:
        raise ValueError('Unexpected step_number={}'.format(args.step_number))

    return output


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: display arrangement of EMIR CSU bars',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # positional arguments
    parser.add_argument("filename",
                        help="TXT file with list of ABBA FITS files",
                        type=argparse.FileType('rt'))
    parser.add_argument("step_number",
                        help=textwrap.dedent("""\
                        0: preliminary rectwv_coeff.json
                        1: refined rectwv_coeff.json
                        2: ABBA reduction"""),
                        type=int, choices=[0, 1, 2])
    parser.add_argument("outfile",
                        help="Output YAML file name",
                        type=lambda x: arg_file_is_new(parser, x))

    # optional arguments
    parser.add_argument("--pattern",
                        help="Observation pattern",
                        default='ABBA',
                        choices=['A', 'AB', 'ABBA'])
    parser.add_argument("--repeat",
                        help="Repetitions of each position",
                        default=1, type=int)
    parser.add_argument("--npreliminary",
                        help="number of images to be combined to compute "
                             "preliminary rectwv_coeff.json",
                        type=int, default=1)
    parser.add_argument("--refine_wavecalib_mode",
                        help=textwrap.dedent("""\
                        0: no refinement
                        1: global offset to all the slitlets (ARC lines)
                        2: individual offset to each slitlet (ARC lines)
                        11: global offset to all the slitlets (OH lines)
                        12: individual offset to each slitlet (OH lines)"""),
                        type=int, choices=[0, 1, 2, 11, 12])
    parser.add_argument("--minimum_slitlet_width_mm",
                        type=float)
    parser.add_argument("--maximum_slitlet_width_mm",
                        type=float)
    parser.add_argument("--global_integer_offset_x_pix",
                        type=int)
    parser.add_argument("--global_integer_offset_y_pix",
                        type=int)
    parser.add_argument("--obsid_prefix",
                        type=str)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")
    args = parser.parse_args(args)

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # read TXT file
    with args.filename as f:
        file_content = f.read().splitlines()
    list_fileinfo = []
    for line in file_content:
        if len(line) > 0:
            if line[0] not in ['#', '@']:
                tmplist = line.split()
                tmpfile = tmplist[0]
                if len(tmplist) > 1:
                    tmpinfo = tmplist[1:]
                else:
                    tmpinfo = None
                list_fileinfo.append(FileInfo(tmpfile, tmpinfo))

    # check consistency of pattern, repeat and number of images
    nimages = len(list_fileinfo)
    pattern_sequence = ''
    for i in range(len(args.pattern)):
        pattern_sequence += args.pattern[i] * args.repeat
    if nimages % len(pattern_sequence) != 0:
        raise ValueError('Unexpected number of images')
    nsequences = nimages // len(pattern_sequence)
    full_set = pattern_sequence * nsequences

    print('Expected sequence pattern: {}'.format(pattern_sequence))
    print('Number of sequences......: {}'.format(nsequences))
    print('Full set of images.......: {}'.format(full_set))

    if args.step_number == 0:
        output = generate_yaml_content(0, args, list_fileinfo)
    elif args.step_number == 1:
        output = generate_yaml_content(1, args, list_fileinfo)
    elif args.step_number == 2:
        output = generate_yaml_content(2, args, list_fileinfo)
    else:
        raise ValueError('Unexpected step_number={}'.format(args.step_number))

    # generate YAML file
    with args.outfile as f:
        f.write(output)
    print('--> File {} generated!'.format(args.outfile.name))


if __name__ == "__main__":
    main()
