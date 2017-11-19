from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
import sys

from numina.array.bpm import process_bpm_median
import numina.array._combine

from arg_file_is_new import arg_file_is_new

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file name",
                        type=argparse.FileType('r'))
    parser.add_argument("--bpm", required=True,
                        help="Bad pixel mask",
                        type=argparse.FileType('r'))
    parser.add_argument("--outfile", required=True,
                        help="Output FITS file name",
                        type=lambda x: arg_file_is_new(parser, x))

    # optional arguments
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting & debugging options"
                             " (default=12)",
                        default=12, type=int,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args()

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # read input FITS file
    hdulist_image = fits.open(args.fitsfile, mode='readonly')
    image_header = hdulist_image[0].header
    image2d = hdulist_image[0].data

    # protections
    naxis1 = image_header['naxis1']
    naxis2 = image_header['naxis2']

    if image2d.shape != (naxis2, naxis1):
        raise ValueError("Unexpected error with NAXIS1, NAXIS2")

    if image2d.shape != (EMIR_NAXIS2, EMIR_NAXIS1):
        raise ValueError("NAXIS1, NAXIS2 unexpected for EMIR detector")

    # read bad pixel mask
    hdulist_bpm = fits.open(args.bpm)
    image_bpm = hdulist_bpm[0].data
    hdulist_bpm.close()

    # apply bad pixel mask
    print('>>', image2d.mean(), image_bpm.mean())
    hdulist_image[0].data = process_bpm_median(image2d, image_bpm)
    print('>>', hdulist_image[0].data.mean(), image_bpm.mean())

    # save output FITS file
    hdulist_image.writeto(args.outfile)

    # close original image
    hdulist_image.close()


if __name__ == "__main__":

    main()
