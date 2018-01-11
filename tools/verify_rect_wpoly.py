from __future__ import division
from __future__ import print_function

import argparse
import astropy.io.fits as fits
import decimal
import json
import numpy as np
import os
from scipy import ndimage
import sys
from uuid import uuid4

from numina.array.display.matplotlib_qt import plt
from numina.array.wavecalib.__main__ import read_wv_master_file
from numina.array.wavecalib.check_wlcalib import check_wlcalib_sp
from numina.array.wavecalib.check_wlcalib import update_poly_wlcalib

from arg_file_is_new import arg_file_is_new
from rect_wpoly_for_mos import islitlet_progress

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NAXIS1


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(prog='verify_rect_wpoly')

    # positional parameters
    parser.add_argument("fitsfile",
                        help="Rectified and wavelength calibrated FITS image",
                        type=argparse.FileType('r'))
    parser.add_argument("--coef_rect_wpoly", required=True,
                        help="Initial JSON file with rectification and "
                             "wavelength calibration coefficients",
                        type=argparse.FileType('r'))
    parser.add_argument("--wv_master_file", required=True,
                        help="TXT file containing wavelengths",
                        type=argparse.FileType('r'))
    parser.add_argument("--verified_coef_rect_wpoly", required=True,
                        help="Output JSON file with improved wavelength "
                             "calibration",
                        type=lambda x: arg_file_is_new(parser, x))

    # optional arguments
    parser.add_argument("--interactive",
                        help="Ask the user for confirmation before updating "
                             "the wavelength calibration polynomial",
                        action="store_true")
    parser.add_argument("--min_nlines",
                        help="Minimum number of lines that must be identified "
                             "(default=0: no minimum is required)",
                        default=0, type=int)
    parser.add_argument("--threshold",
                        help="Minimum signal in the line peaks (default=0)",
                        default=0, type=float)
    parser.add_argument("--sigma_gaussian_filtering",
                        help="Sigma of the gaussian filter to be applied to "
                             "the spectrum in order to avoid problems with "
                             "saturated lines",
                        default=0, type=float)
    parser.add_argument("--nwinwidth_initial",
                        help="Width of the window (pixels) where each peak "
                             "must be initially found (default 7)",
                        default=7, type=int)
    parser.add_argument("--nwinwidth_refined",
                        help="Width of the window (pixels) where each peak "
                             "must be refined (default 5)",
                        default=5, type=int)
    parser.add_argument("--ntimes_match_wv",
                        help="Times CDELT1 to match measured and expected "
                             "wavelengths (default 2)",
                        default=2, type=float)
    parser.add_argument("--poldeg_residuals",
                        help="Polynomial degree for fit to residuals "
                             "(default 1)",
                        default=1, type=int)
    parser.add_argument("--times_sigma_reject",
                        help="Times the standard deviation to reject points "
                             "iteratively in the fit to residuals ("
                             "default=5)",
                        default=5, type=float)
    parser.add_argument("--use_r",
                        help="Perform additional statistical analysis with R",
                        action="store_true")
    parser.add_argument("--geometry",
                        help="tuple x,y,dx,dy (default 0,0,640,480)",
                        default="0,0,640,480")
    parser.add_argument("--debugplot",
                        help="integer indicating plotting/debugging" +
                             " (default=-11)",
                        type=int, default=-11,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args(args=args)

    if args.echo:
        print('\033[1m\033[31m% ' + ' '.join(sys.argv) + '\033[0m\n')

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

    # read calibration structure from JSON file
    rect_wpoly_dict = json.loads(open(args.coef_rect_wpoly.name).read())
    islitlet_min = rect_wpoly_dict['tags']['islitlet_min']
    islitlet_max = rect_wpoly_dict['tags']['islitlet_max']

    # read FITS image and its corresponding header
    hdulist = fits.open(args.fitsfile)
    header = hdulist[0].header
    image2d_rectified_wv = hdulist[0].data
    hdulist.close()
    crpix1_enlarged = header['crpix1']
    crval1_enlarged = header['crval1']
    cdelt1_enlarged = header['cdelt1']

    # read master arc line wavelengths (whole data set)
    wv_master_all = read_wv_master_file(
        wv_master_file=args.wv_master_file.name,
        lines='all',
        debugplot=args.debugplot
    )

    basefilename = os.path.basename(args.fitsfile.name)

    # main loop
    for islitlet in range(islitlet_min, islitlet_max + 1):
        if args.debugplot == 0 and not args.interactive:
            islitlet_progress(islitlet, islitlet_max)
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        if abs(args.debugplot) >= 10:
            print('=' * 79)
            print(">>> " + cslitlet)

        # scan region spanned by current slitlet
        sltmin = header['sltmin' + str(islitlet).zfill(2)]
        sltmax = header['sltmax' + str(islitlet).zfill(2)]

        # median spectrum
        spmedian = np.median(
            image2d_rectified_wv[sltmin:(sltmax + 1)],
            axis=0
        )

        # gaussian filtering when requested (to avoid line saturation)
        if args.sigma_gaussian_filtering > 0:
            spmedian = ndimage.filters.gaussian_filter(
                spmedian,
                sigma=args.sigma_gaussian_filtering
            )

        # check wavelength calibration
        polyres, ysummary, xyrfit = check_wlcalib_sp(
            sp=spmedian,
            crpix1=crpix1_enlarged,
            crval1=crval1_enlarged,
            cdelt1=cdelt1_enlarged,
            wv_master=wv_master_all,
            threshold=args.threshold,
            nwinwidth_initial=args.nwinwidth_initial,
            nwinwidth_refined=args.nwinwidth_refined,
            ntimes_match_wv=args.ntimes_match_wv,
            poldeg_residuals=args.poldeg_residuals,
            times_sigma_reject=args.times_sigma_reject,
            use_r=args.use_r,
            title=basefilename + ' [slitlet #' + str(islitlet).zfill(2) + ']',
            full=True,
            geometry=geometry,
            debugplot=abs(args.debugplot))

        npoints_total = len(xyrfit[0])
        npoints_removed = sum(xyrfit[2])
        npoints_used = npoints_total - npoints_removed
        if abs(args.debugplot) >= 10:
            print('>>> Npoints (total / used / removed)..:',
                  npoints_total, npoints_used, npoints_removed)
        if npoints_used > args.min_nlines:
            if args.interactive:
                # include correction to initial wavelength calibration after
                # confirmation by the user
                print("Update wavelength calibration polynomial with "
                      "last fit:")
                loop = True
                while(loop):
                    print("(y)es, (n)o, (u)ser defined correction polynomial "
                          "[y]? ", end='')
                    coption = sys.stdin.readline().strip()
                    coption = coption.lower()
                    if coption == '':
                        coption = 'y'
                    if coption in 'ynu':
                        loop = False
            else:
                coption = 'y'
            if coption == 'y' or coption == 'u':
                wpoly_coeff = rect_wpoly_dict['contents'][cslitlet]['wpoly_coeff']
                if coption == 'y':
                    coeff_residuals = polyres.coef
                elif coption == 'u':
                    print("Polynomial degree? ", end='')
                    degcorr = int(sys.stdin.readline().strip())
                    coeff_residuals = np.zeros(degcorr + 1)
                    for idum in range(degcorr + 1):
                        print("Coefficient #" + str(idum) + "? ", end='')
                        coeff_residuals[idum] = \
                            float(sys.stdin.readline().strip())
                else:
                    raise ValueError('Unexpected coption=' + str(coption))

                wpoly_coeff_updated = update_poly_wlcalib(
                    coeff_ini=wpoly_coeff,
                    coeff_residuals=coeff_residuals,
                    xyrfit=xyrfit,
                    naxis1=EMIR_NAXIS1,
                    debugplot=0
                ).tolist()
                rect_wpoly_dict['contents'][cslitlet]['wpoly_coeff'] = \
                    wpoly_coeff_updated
                if abs(args.debugplot) >= 10:
                    for idum, fdum in \
                            enumerate(zip(wpoly_coeff, wpoly_coeff_updated)):
                        print(">>> coef#" + str(idum) + ':  ', end='')
                        print("%+.8E  -->  %+.8E" % (decimal.Decimal(fdum[0]),
                                                     decimal.Decimal(fdum[1])))
        else:
            if abs(args.debugplot) >= 10:
                print('>>> SKIPPING slitlet: number of lines < min_nlines')
        if args.debugplot < 0:
            plt.close()

    # update uuid for verified JSON structure
    rect_wpoly_dict['uuid'] = str(uuid4())

    # save verified JSON file
    with open(args.verified_coef_rect_wpoly.name, 'w') as fstream:
        json.dump(rect_wpoly_dict, fstream, indent=2, sort_keys=True)
        print('>>> Saving file ' + args.verified_coef_rect_wpoly.name)


if __name__ == "__main__":

    main()
