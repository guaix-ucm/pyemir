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
def generate_yaml_content(args, list_fileinfo, enabled=True):

    step = args.step

    nimages = len(list_fileinfo)

    output = ''
    if step in [0, 1]:
        # initial reduction with STARE_IMAGE
        for i in range(nimages):
            if i % args.repeat == 0:
                if i != 0:
                    output += '---\n'
                idlabel = list_fileinfo[i].filename[:10]
                output += f'id: _{idlabel}'
                output += '\n'
                output += 'instrument: EMIR\n'
                output += 'mode: STARE_IMAGE\n'
                output += 'frames:\n'
            output += f' - {list_fileinfo[i].filename}\n'
            if (i + 1) % args.repeat == 0:
                output += 'requirements:\n'
                output += f'  reprojection_method: {args.reprojection}\n'
                if step == 0:
                    output += 'enabled: True\n'
                elif step == 1:
                    output += 'enabled: False\n'
                else:
                    raise ValueError(f'Unexpected step={step}')

    if step == 1:
        # combination with FULL_DITHERED_IMAGE
        output += '---\n'
        if args.obsid_combined is None:
            output += 'id: _combined' + '\n'
        else:
            output += f'id: _{args.obsid_combined}\n'
        output += 'instrument: EMIR\n'
        output += 'mode: FULL_DITHERED_IMAGE\n'
        output += 'children:\n'
        for i in range(nimages):
            if i % args.repeat == 0:
                idlabel = list_fileinfo[i].filename[:10]
                output += f' - _{idlabel}\n'
        output += 'requirements:\n'
        output += '  iterations: 0\n'
        output += '  sky_images: 0\n'
        output += '  refine_offsets: False\n'
        output += 'enabled: True\n'

    return output


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: generate observation result YAML file',
    )

    # positional or required arguments
    parser.add_argument("filename",
                        help="TXT file with list of ABBA FITS files",
                        type=argparse.FileType('rt'))
    parser.add_argument("--step", required=True,
                        help=textwrap.dedent("""\
                        0: preliminary STARE_IMAGE
                        1: combination with FULL_DITHERED_IMAGE"""),
                        type=int, choices=[0, 1])
    parser.add_argument("--outfile", required=True,
                        help="Output YAML file name",
                        type=lambda x: arg_file_is_new(parser, x))

    # optional arguments
    parser.add_argument("--reprojection",
                        help="Reprojection method",
                        default="interp", type=str,
                        choices=["interp", "adaptive", "exact", "none"])
    parser.add_argument("--repeat",
                        help="Repetitions at each position",
                        default=1, type=int)
    parser.add_argument("--obsid_combined",
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
    if nimages % args.repeat != 0:
        raise ValueError('Unexpected number of images')

    output = generate_yaml_content(args, list_fileinfo)

    # generate YAML file
    with args.outfile as f:
        f.write(output)
    print('--> File {} generated!'.format(args.outfile.name))


if __name__ == "__main__":
    main()
