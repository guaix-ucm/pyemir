#
# Copyright 2018 Universidad Complutense de Madrid
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

from emirdrp.processing.wavecal.rectwv import rectwv_coeff_from_arc_image
from emirdrp.products import RefinedBoundaryModelParam

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from numina.tools.arg_file_is_new import arg_file_is_new


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file with longslit data",
                        type=argparse.FileType('rb'))
    parser.add_argument("--bound_param", required=True,
                        help="Input JSON with fitted boundary parameters",
                        type=argparse.FileType('rt'))
    parser.add_argument("--order_fmap", required=True,
                        help="Order of the 2D rectification transformation "
                             "(default=2)",
                        default=2, type=int)
    parser.add_argument("--wv_master_file", required=True,
                        help="TXT file containing wavelengths")
    parser.add_argument("--poldeg_initial", required=True,
                        help="Polynomial degree for initial calibration",
                        type=int)
    parser.add_argument("--poldeg_refined", required=True,
                        help="Polynomial degree for refined calibration "
                             "(0=do not refine)",
                        type=int)
    parser.add_argument("--out_json", required=True,
                        help="Output JSON file with results",
                        type=lambda x: arg_file_is_new(parser, x))

    # optional arguments
    parser.add_argument("--interactive",
                        help="Ask the user for confirmation before updating "
                             "the wavelength calibration polynomial",
                        action="store_true")
    parser.add_argument("--ymargin_bb",
                        help="Number of pixels above and below frontiers to "
                             "determine the vertical bounding box of each "
                             "undistorted slitlet (default=2)",
                        type=int, default=2)
    parser.add_argument("--remove_sp_background",
                        help="Remove background spectrum prior to arc line "
                             "detection",
                        action="store_true")
    parser.add_argument("--times_sigma_threshold",
                        help="Times sigma above threshold to detect unknown"
                             " arc lines (default=10)",
                        type=float, default=10)
    parser.add_argument("--margin_npix",
                        help="Number of pixels before and after expected "
                             "wavelength calibrated spectrum to trim the "
                             "wv_master table in the wavelength direction "
                             "(default=50)",
                        type=int, default=50)
    parser.add_argument("--nbrightlines",
                        help="tuple with number of brightlines to "
                             "be employed in the initial wavelength "
                             "calibration (e.g. \"10,5,4\")")
    parser.add_argument("--threshold_wv",
                        help="Minimum signal in the line peaks (default=0)",
                        default=0, type=float)
    parser.add_argument("--sigma_gaussian_filtering",
                        help="Sigma of the gaussian filter to be applied to "
                             "the spectrum in order to avoid problems with "
                             "saturated lines in the wavelength calibration "
                             "process",
                        default=0, type=float)
    parser.add_argument("--out_55sp",
                        help="FITS file containing the set of averaged "
                             "spectra employed to derive the wavelength "
                             "calibration",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))
    parser.add_argument("--ylogscale",
                        help="Display spectrum signal in logarithmic units",
                        action="store_true")
    parser.add_argument("--geometry",
                        help="tuple x,y,dx,dy (default 0,0,640,480)",
                        default="0,0,640,480")
    parser.add_argument("--pdffile",
                        help="output PDF file name",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))
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

    # read pdffile
    if args.pdffile is not None:
        if args.interactive:
            raise ValueError('--interactive is not compatible with --pdffile')
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages(args.pdffile.name)
    else:
        pdf = None

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

    # generate RefinedBoundaryModelParam object
    bound_param = RefinedBoundaryModelParam._datatype_load(
        args.bound_param.name)

    # generate HDUList object
    hdulist = fits.open(args.fitsfile)

    # generate lines_catalog
    lines_catalog = np.genfromtxt(args.wv_master_file)

    rectwv_coeff, reduced_55sp = rectwv_coeff_from_arc_image(
        hdulist,
        bound_param,
        lines_catalog,
        args_nbrightlines=args.nbrightlines,
        args_ymargin_bb=args.ymargin_bb,
        args_remove_sp_background=args.remove_sp_background,
        args_times_sigma_threshold=args.times_sigma_threshold,
        args_order_fmap=args.order_fmap,
        args_sigma_gaussian_filtering=args.sigma_gaussian_filtering,
        args_margin_npix=args.margin_npix,
        args_poldeg_initial=args.poldeg_initial,
        args_poldeg_refined=args.poldeg_refined,
        args_interactive=args.interactive,
        args_threshold_wv=args.threshold_wv,
        args_ylogscale=args.ylogscale,
        args_pdf=pdf,
        args_geometry=geometry,
        debugplot=0
    )

    # Save image with collapsed spectra employed to determine the
    # wavelength calibration
    if args.out_55sp is not None:
        reduced_55sp.writeto(args.out_55sp, overwrite=True)

    # save RectWaveCoeff object into JSON file
    rectwv_coeff.writeto(args.out_json.name)
    print('>>> Saving file ' + args.out_json.name)
    # debugging __getstate__ and __setstate__
    # check_setstate_getstate(rectwv_coeff, args.out_json.name)

    if pdf is not None:
        pdf.close()


if __name__ == "__main__":
    main()
