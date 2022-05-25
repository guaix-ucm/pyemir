#
# Copyright 2008-2020 Universidad Complutense de Madrid
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
import numpy as np
import sys

from numina.array.wavecalib.fix_pix_borders import find_pix_borders
from numina.frame.utils import copy_img
from numina.tools.arg_file_is_new import arg_file_is_new

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def merge2images(hdu1, hdu2, debugplot):
    """Merge 2 EMIR images, averaging the common region.

    Parameters
    ----------
    hdu1 : HDUList object
        Input image #1.
    hdu2 : HDUList object
        Input image #2.
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed. The valid codes are defined in
        numina.array.display.pause_debugplot.

    Returns
    -------
    image_merged : HDUList object
        Output image.

    """

    # check image dimensions
    image_header = hdu1[0].header
    image_merged = copy_img(hdu1)
    image_header_ = hdu2[0].header
    #
    naxis = image_header['naxis']
    naxis_ = image_header_['naxis']
    if naxis != naxis_:
        raise ValueError('Incompatible NAXIS values: {}, {}'.format(
            naxis, naxis_
        ))
    naxis1 = image_header['naxis1']
    naxis1_ = image_header_['naxis1']
    if naxis1 != naxis1_:
        raise ValueError('Incompatible NAXIS1 values: {}, {}'.format(
            naxis1, naxis1_
        ))
    if naxis == 1:
        naxis2 = 1
    elif naxis == 2:
        naxis2 = image_header['naxis2']
        naxis2_ = image_header_['naxis2']
        if naxis2 != naxis2_:
            raise ValueError('Incompatible NAXIS2 values: {}, {}'.format(
                naxis2, naxis2_
            ))
    else:
        raise ValueError('Unexpected NAXIS value: {}'.format(naxis))

    # initialize output array
    image2d_merged = np.zeros((naxis2, naxis1))

    # main loop
    data1 = hdu1[0].data
    data2 = hdu2[0].data
    for i in range(naxis2):
        sp1 = data1[i, :]
        sp2 = data2[i, :]
        jmin1, jmax1 = find_pix_borders(sp1, sought_value=0)
        jmin2, jmax2 = find_pix_borders(sp2, sought_value=0)
        image2d_merged[i, :] = sp1 + sp2
        useful1 = (jmin1 != -1) and (jmax1 != naxis1)
        useful2 = (jmin2 != -1) and (jmax2 != naxis1)
        if useful1 and useful2:
            jmineff = max(jmin1, jmin2)
            jmaxeff = min(jmax1, jmax2)
            image2d_merged[i, jmineff:(jmaxeff+1)] /= 2

    # return result
    image_merged[0].data = image2d_merged.astype(np.float32)
    return image_merged


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: merge 2 EMIR images averaging the common'
                    ' region'
    )

    # positional arguments
    parser.add_argument("infile1",
                        help="Input FITS file name #1",
                        type=argparse.FileType('rb'))
    parser.add_argument("infile2",
                        help="Input FITS file name #2",
                        type=argparse.FileType('rb'))
    parser.add_argument("outfile",
                        help="Output FITS file name",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))

    # optional arguments
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting/debugging" +
                             " (default=0)",
                        default=0, type=int,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args(args=args)

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # read input FITS files
    hdulist1 = fits.open(args.infile1)
    hdulist2 = fits.open(args.infile2)

    image_merged = merge2images(
        hdulist1,
        hdulist2,
        debugplot=args.debugplot
    )

    # save result
    image_merged.writeto(args.outfile, overwrite=True)


if __name__ == "__main__":

    main()