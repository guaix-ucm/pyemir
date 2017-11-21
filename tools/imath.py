#
# Copyright 2015-2017 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# Numina is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Numina is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Numina.  If not, see <http://www.gnu.org/licenses/>.
#

"""Auxiliary script to perform basic image arithmetic."""

from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
import numpy as np
import sys

from numina.array.display.ximshow import ximshow_file

from arg_file_is_new import arg_file_is_new


def compute_operation(file1, file2, operation, output, display,
                      args_z1z2, args_bbox, args_keystitle, args_geometry):
    """Compute output = file1 operation file2.

    Parameters
    ----------
    file1 : string
        First FITS file.
    file2 : string
        Second FITS file.
    operation : string
        Mathematical operation.
    output : string
        Output FITS file.
    display : string
        Character string indication whether the images are displayed.
        Valid values are 'all', 'result' and 'none' (default).
    args_z1z2 : string or None
        String providing the image cuts tuple: z1, z2, minmax or None.
    args_bbox : string or None
        String providing the bounding box tuple: nc1, nc2, ns1, ns2.
    args_keystitle : string or None
        Tuple of FITS keywords.format: key1,key2,...,keyn.format
    args_geometry : string or None
        Tuple x,y,dx,dy to define the Qt4 backend geometry.

    """

    # read first FITS file
    with fits.open(file1) as hdulist:
        image_header1 = hdulist[0].header
        image1 = hdulist[0].data.astype(np.float)
    naxis1 = image_header1['naxis1']
    naxis2 = image_header1['naxis2']

    # if required, display file1
    if display == 'all':
        ximshow_file(file1,
                     args_z1z2=args_z1z2, args_bbox=args_bbox,
                     args_keystitle=args_keystitle,
                     args_geometry=args_geometry,
                     debugplot=12)

    # read second FITS file
    with fits.open(file2) as hdulist:
        image_header2 = hdulist[0].header
        image2 = hdulist[0].data.astype(np.float)
    naxis1_ = image_header2['naxis1']
    naxis2_ = image_header2['naxis2']

    # if required, display file2
    if display == 'all':
        ximshow_file(file2,
                     args_z1z2=args_z1z2, args_bbox=args_bbox,
                     args_keystitle=args_keystitle,
                     args_geometry=args_geometry,
                     debugplot=12)

    # check dimensions
    if naxis1 != naxis1_:
        raise ValueError("NAXIS1 values are different.")
    if naxis2 != naxis2_:
        raise ValueError("NAXIS2 values are different.")

    # compute operation
    if operation == "+":
        solution = image1 + image2
    elif operation == "-":
        solution = image1 - image2
    elif operation == "*":
        solution = image1 * image2
    elif operation == "/":
        solution = image1 / image2
    else:
        raise ValueError("Unexpected operation=" + str(operation))

    # save output file
    hdu = fits.PrimaryHDU(solution.astype(np.float), image_header1)
    hdu.writeto(output, overwrite=True)

    # if required, display result
    if display in ['all', 'result']:
        ximshow_file(output,
                     args_z1z2=args_z1z2, args_bbox=args_bbox,
                     args_keystitle=args_keystitle,
                     args_geometry=args_geometry,
                     debugplot=12)


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()
    # positional parameters
    parser.add_argument("file1",
                        help="First FITS image",
                        type=argparse.FileType('r'))
    parser.add_argument("operation",
                        help="Mathematical operation",
                        type=str,
                        choices=['+', '-', '*', '/'])
    parser.add_argument("file2",
                        help="Second FITS image",
                        type=argparse.FileType('r'))
    # optional arguments
    parser.add_argument("output",
                        help="Output FITS image",
                        type=lambda x: arg_file_is_new(parser, x))
    parser.add_argument("--display",
                        help="Display images: all, result, none (default)",
                        default="none",
                        type=str,
                        choices=['all', 'result', 'none'])
    parser.add_argument("--z1z2",
                        help="tuple z1,z2, minmax or None (use zscale)")
    parser.add_argument("--bbox",
                        help="bounding box tuple: nc1,nc2,ns1,ns2")
    parser.add_argument("--keystitle",
                        help="tuple of FITS keywords.format: " +
                             "key1,key2,...keyn.'format'")
    parser.add_argument("--geometry",
                        help="tuple x,y,dx,dy")
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args(args=args)

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # compute operation
    compute_operation(args.file1, args.file2,
                      args.operation, args.output,
                      args.display,
                      args.z1z2, args.bbox, args.keystitle, args.geometry)


if __name__ == "__main__":

    main()
