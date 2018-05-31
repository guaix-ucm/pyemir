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
import numpy as np
import sys

from numina.tools.arg_file_is_new import arg_file_is_new

from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_NPIXPERSLIT_RECTIFIED


def median_slitlets_rectified(input_image, sp55=False):
    """Compute median spectrum for each slitlet

    Parameters
    ----------
    input_image : HDUList object
        Input 2D image.
    sp55 : bool
        If True, the returned array contains only EMIR_NBARS
        spectra, one for each slitlet. If False, the returned
        image has the same dimensions as the input image, with the
        median spectrum of each slitlet spanning all the spectra
        of the corresponding slitlet.

    Returns
    -------
    image_median : HDUList object
        Output 2d image where the spectra of each slitlet has been
        replaced by the median spectrum of the same slitlet.

    """

    image_header = input_image[0].header
    image2d = input_image[0].data

    # check image dimensions
    naxis2_expected = EMIR_NBARS * EMIR_NPIXPERSLIT_RECTIFIED

    naxis2, naxis1 = image2d.shape
    if naxis2 != naxis2_expected:
        raise ValueError("NAXIS2={0} should be {1}".format(
            naxis2, naxis2_expected
        ))

    # check that the FITS file has been obtained with EMIR
    instrument = image_header['instrume']
    if instrument != 'EMIR':
        raise ValueError("INSTRUME keyword is not 'EMIR'!")

    # initialize output image
    if sp55:
        image2d_median = np.zeros((EMIR_NBARS, naxis1))
    else:
        image2d_median = np.zeros((naxis2, naxis1))

    # main loop
    for i in range(EMIR_NBARS):
        ns1 = i * EMIR_NPIXPERSLIT_RECTIFIED + 1
        ns2 = ns1 + EMIR_NPIXPERSLIT_RECTIFIED - 1
        sp_median = np.median(image2d[(ns1-1):ns2, :], axis=0)

        if sp55:
            image2d_median[i] = np.copy(sp_median)
        else:
            image2d_median[(ns1-1):ns2, :] = np.tile(
                sp_median, (EMIR_NPIXPERSLIT_RECTIFIED, 1)
            )

    image_median = fits.PrimaryHDU(data=image2d_median, header=image_header)

    return image_median


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: overplot boundary model over FITS image'
    )

    # positional arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file name",
                        type=argparse.FileType('rb'))
    parser.add_argument("outfile",
                        help="Output FITS file name",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))

    # optional arguments
    parser.add_argument("--only55",
                        help="Do not expand the median spectrum of each "
                             "slitlet",
                        action="store_true")
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args()

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # read input FITS file
    hdulist = fits.open(args.fitsfile)

    image_median = median_slitlets_rectified(
        hdulist,
        sp55=args.only55
    )

    # save result
    image_median.writeto(args.outfile, overwrite=True)


if __name__ == "__main__":

    main()
