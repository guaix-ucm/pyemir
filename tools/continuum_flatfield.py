from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
import json
import numpy as np
import sys

from numina.array.wavecalib.resample import resample_image2d_flux

from apply_rect_wpoly import Slitlet2D
from arg_file_is_new import arg_file_is_new
from emirdrp.instrument.dtu_configuration import DtuConfiguration
from nscan_minmax_frontiers import nscan_minmax_frontiers
from rect_wpoly_for_mos import islitlet_progress
from save_ndarray_to_fits import save_ndarray_to_fits

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(prog='apply_rect_wpoly')

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file (flat ON-OFF)",
                        type=argparse.FileType('r'))
    parser.add_argument("--coef_rect_wpoly", required=True,
                        help="Input JSON file with rectification and "
                             "wavelength calibration coefficients",
                        type=argparse.FileType('r'))
    parser.add_argument("--outfile", required=True,
                        help="Output FITS file",
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

    # read calibration structure from JSON file
    rect_wpoly_dict = json.loads(open(args.coef_rect_wpoly.name).read())

    # read FITS image and its corresponding header
    hdulist = fits.open(args.fitsfile)
    header = hdulist[0].header
    image2d = hdulist[0].data
    hdulist.close()

    # protections
    naxis2, naxis1 = image2d.shape
    if naxis1 != header['naxis1'] or naxis2 != header['naxis2']:
        print('>>> NAXIS1:', naxis1)
        print('>>> NAXIS2:', naxis2)
        raise ValueError('Something is wrong with NAXIS1 and/or NAXIS2')
    if abs(args.debugplot) >= 10:
        print('>>> NAXIS1:', naxis1)
        print('>>> NAXIS2:', naxis2)

    # check that the input FITS file grism and filter match
    filter_name = header['filter']
    if filter_name != rect_wpoly_dict['tags']['filter']:
        raise ValueError("Filter name does not match!")
    grism_name = header['grism']
    if grism_name != rect_wpoly_dict['tags']['grism']:
        raise ValueError("Filter name does not match!")
    if abs(args.debugplot) >= 10:
        print('>>> grism.......:', grism_name)
        print('>>> filter......:', filter_name)

    # check that the DTU configurations are compatible
    dtu_conf_fitsfile = DtuConfiguration()
    dtu_conf_fitsfile.define_from_fits(args.fitsfile)
    dtu_conf_jsonfile = DtuConfiguration()
    dtu_conf_jsonfile.define_from_dictionary(
        rect_wpoly_dict['dtu_configuration'])
    if dtu_conf_fitsfile != dtu_conf_jsonfile:
        print('DTU configuration (FITS file):\n\t', dtu_conf_fitsfile)
        print('DTU configuration (JSON file):\n\t', dtu_conf_jsonfile)
        raise ValueError('Incompatible DTU configurations')

    # read islitlet_min and islitlet_max from input JSON file
    islitlet_min = rect_wpoly_dict['tags']['islitlet_min']
    islitlet_max = rect_wpoly_dict['tags']['islitlet_max']
    if abs(args.debugplot) >= 10:
        print('>>> islitlet_min:', islitlet_min)
        print('>>> islitlet_max:', islitlet_max)

    # ---

    # initialize rectified image
    image2d_flatfielded = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

    # main loop
    for islitlet in range(islitlet_min, islitlet_max + 1):
        if args.debugplot == 0:
            islitlet_progress(islitlet, islitlet_max)

        # define Slitlet2D object
        slt = Slitlet2D(islitlet=islitlet,
                        megadict=rect_wpoly_dict,
                        ymargin=2,
                        debugplot=args.debugplot)

        if abs(args.debugplot) >= 10:
            print(slt)

        # extract (distorted) slitlet from the initial image
        slitlet2d = slt.extract_slitlet2d(image2d)

        # rectify slitlet
        slitlet2d_rect = slt.rectify(slitlet2d, resampling=1)

        # ToDo: compute median spectrum, normalize, and unrectify image

    save_ndarray_to_fits(
        array=image2d_flatfielded,
        file_name=args.outfile,
        main_header=header,
        overwrite=True
    )
    print('>>> Saving file ' + args.outfile.name)


if __name__ == "__main__":
    main()
