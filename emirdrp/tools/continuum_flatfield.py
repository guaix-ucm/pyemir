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
from scipy import ndimage
import sys

from numina.array.display.ximplotxy import ximplotxy
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.wavecalib.apply_integer_offsets import apply_integer_offsets
from numina.tools.arg_file_is_new import arg_file_is_new

from emirdrp.processing.wavecal.slitlet2d import Slitlet2D
from emirdrp.instrument.dtu_configuration import DtuConfiguration
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers
from emirdrp.tools.rect_wpoly_for_mos import islitlet_progress
from emirdrp.tools.save_ndarray_to_fits import save_ndarray_to_fits
from emirdrp.products import RectWaveCoeff

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: compute pixel-to-pixel flatfield'
    )

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file (flat ON-OFF)",
                        type=argparse.FileType('rb'))
    parser.add_argument("--rectwv_coeff", required=True,
                        help="Input JSON file with rectification and "
                             "wavelength calibration coefficients",
                        type=argparse.FileType('rt'))
    parser.add_argument("--minimum_fraction", required=True,
                        help="Minimum allowed flatfielding value",
                        type=float, default=0.01)
    parser.add_argument("--minimum_value_in_output",
                        help="Minimum value allowed in output file: pixels "
                             "below this value are set to 1.0 (default=0.01)",
                        type=float, default=0.01)
    parser.add_argument("--nwindow_median", required=True,
                        help="Window size to smooth median spectrum in the "
                             "spectral direction",
                        type=int)
    parser.add_argument("--outfile", required=True,
                        help="Output FITS file",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))

    # optional arguments
    parser.add_argument("--delta_global_integer_offset_x_pix",
                        help="Delta global integer offset in the X direction "
                             "(default=0)",
                        default=0, type=int)
    parser.add_argument("--delta_global_integer_offset_y_pix",
                        help="Delta global integer offset in the Y direction "
                             "(default=0)",
                        default=0, type=int)
    parser.add_argument("--resampling",
                        help="Resampling method: 1 -> nearest neighbor, "
                             "2 -> linear interpolation (default)",
                        default=2, type=int,
                        choices=(1, 2))
    parser.add_argument("--ignore_DTUconf",
                        help="Ignore DTU configurations differences between "
                             "model and input image",
                        action="store_true")
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
    rectwv_coeff = RectWaveCoeff._datatype_load(args.rectwv_coeff.name)

    # modify (when requested) global offsets
    rectwv_coeff.global_integer_offset_x_pix += \
        args.delta_global_integer_offset_x_pix
    rectwv_coeff.global_integer_offset_y_pix += \
        args.delta_global_integer_offset_y_pix

    # read FITS image and its corresponding header
    hdulist = fits.open(args.fitsfile)
    header = hdulist[0].header
    image2d = hdulist[0].data
    hdulist.close()

    # apply global offsets
    image2d = apply_integer_offsets(
        image2d=image2d,
        offx=rectwv_coeff.global_integer_offset_x_pix,
        offy=rectwv_coeff.global_integer_offset_y_pix
    )

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
    if filter_name != rectwv_coeff.tags['filter']:
        raise ValueError("Filter name does not match!")
    grism_name = header['grism']
    if grism_name != rectwv_coeff.tags['grism']:
        raise ValueError("Filter name does not match!")
    if abs(args.debugplot) >= 10:
        print('>>> grism.......:', grism_name)
        print('>>> filter......:', filter_name)

    # check that the DTU configurations are compatible
    dtu_conf_fitsfile = DtuConfiguration.define_from_fits(args.fitsfile)
    dtu_conf_jsonfile = DtuConfiguration.define_from_dictionary(
        rectwv_coeff.meta_info['dtu_configuration'])
    if dtu_conf_fitsfile != dtu_conf_jsonfile:
        print('DTU configuration (FITS file):\n\t', dtu_conf_fitsfile)
        print('DTU configuration (JSON file):\n\t', dtu_conf_jsonfile)
        if args.ignore_DTUconf:
            print('WARNING: DTU configuration differences found!')
        else:
            raise ValueError('DTU configurations do not match')
    else:
        if abs(args.debugplot) >= 10:
            print('>>> DTU Configuration match!')
            print(dtu_conf_fitsfile)

    # valid slitlet numbers
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in rectwv_coeff.missing_slitlets:
        list_valid_islitlets.remove(idel)
    if abs(args.debugplot) >= 10:
        print('>>> valid slitlet numbers:\n', list_valid_islitlets)

    # ---

    # initialize rectified image
    image2d_flatfielded = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

    # main loop
    for islitlet in list_valid_islitlets:
        if args.debugplot == 0:
            islitlet_progress(islitlet, EMIR_NBARS)

        # define Slitlet2D object
        slt = Slitlet2D(islitlet=islitlet,
                        rectwv_coeff=rectwv_coeff,
                        debugplot=args.debugplot)

        if abs(args.debugplot) >= 10:
            print(slt)

        # extract (distorted) slitlet from the initial image
        slitlet2d = slt.extract_slitlet2d(image2d)

        # rectify slitlet
        slitlet2d_rect = slt.rectify(
            slitlet2d,
            resampling=args.resampling
        )
        naxis2_slitlet2d, naxis1_slitlet2d = slitlet2d_rect.shape

        if naxis1_slitlet2d != EMIR_NAXIS1:
            print('naxis1_slitlet2d: ', naxis1_slitlet2d)
            print('EMIR_NAXIS1.....: ', EMIR_NAXIS1)
            raise ValueError("Unexpected naxis1_slitlet2d")

        # get useful slitlet region (use boundaires instead of frontiers;
        # note that the nscan_minmax_frontiers() works well independently
        # of using frontiers of boundaries as arguments)
        nscan_min, nscan_max = nscan_minmax_frontiers(
            slt.y0_reference_lower,
            slt.y0_reference_upper,
            resize=False
        )
        ii1 = nscan_min - slt.bb_ns1_orig
        ii2 = nscan_max - slt.bb_ns1_orig + 1

        # median spectrum
        sp_collapsed = np.median(slitlet2d_rect[ii1:(ii2 + 1), :], axis=0)

        # smooth median spectrum along the spectral direction
        sp_median = ndimage.median_filter(sp_collapsed, args.nwindow_median,
                                          mode='nearest')
        ymax_spmedian = sp_median.max()
        y_threshold = ymax_spmedian * args.minimum_fraction
        sp_median[np.where(sp_median < y_threshold)] = 0.0

        if abs(args.debugplot) > 10:
            title = 'Slitlet#' + str(islitlet) + '(median spectrum)'
            xdum = np.arange(1, naxis1_slitlet2d + 1)
            ax = ximplotxy(xdum, sp_collapsed,
                           title=title,
                           show=False, **{'label' : 'collapsed spectrum'})
            ax.plot(xdum, sp_median, label='filtered spectrum')
            ax.plot([1, naxis1_slitlet2d], 2*[y_threshold],
                    label='threshold')
            ax.legend()
            ax.set_ylim(-0.05*ymax_spmedian, 1.05*ymax_spmedian)
            pause_debugplot(args.debugplot,
                            pltshow=True, tight_layout=True)

        # generate rectified slitlet region filled with the median spectrum
        slitlet2d_rect_spmedian = np.tile(sp_median, (naxis2_slitlet2d, 1))
        if abs(args.debugplot) > 10:
            slt.ximshow_rectified(slitlet2d_rect_spmedian)

        # unrectified image
        slitlet2d_unrect_spmedian = slt.rectify(
            slitlet2d_rect_spmedian,
            resampling=args.resampling,
            inverse=True
        )

        # normalize initial slitlet image (avoid division by zero)
        slitlet2d_norm = np.zeros_like(slitlet2d)
        for j in range(naxis1_slitlet2d):
            for i in range(naxis2_slitlet2d):
                den = slitlet2d_unrect_spmedian[i, j]
                if den == 0:
                    slitlet2d_norm[i, j] = 1.0
                else:
                    slitlet2d_norm[i, j] = slitlet2d[i, j] / den

        if abs(args.debugplot) > 10:
            slt.ximshow_unrectified(slitlet2d_norm)

        for j in range(EMIR_NAXIS1):
            xchannel = j + 1
            y0_lower = slt.list_frontiers[0](xchannel)
            y0_upper = slt.list_frontiers[1](xchannel)
            n1, n2 = nscan_minmax_frontiers(y0_frontier_lower=y0_lower,
                                            y0_frontier_upper=y0_upper,
                                            resize=True)
            # note that n1 and n2 are scans (ranging from 1 to NAXIS2)
            nn1 = n1 - slt.bb_ns1_orig + 1
            nn2 = n2 - slt.bb_ns1_orig + 1
            image2d_flatfielded[(n1 - 1):n2, j] = \
                slitlet2d_norm[(nn1 - 1):nn2, j]

            # force to 1.0 region around frontiers
            image2d_flatfielded[(n1 - 1):(n1 + 2), j] = 1
            image2d_flatfielded[(n2 - 5):n2, j] = 1
    if args.debugplot == 0:
        print('OK!')

    # set pixels below minimum value to 1.0
    filtered = np.where(image2d_flatfielded < args.minimum_value_in_output)
    image2d_flatfielded[filtered] = 1.0

    # restore global offsets
    image2d_flatfielded = apply_integer_offsets(
        image2d=image2d_flatfielded ,
        offx=-rectwv_coeff.global_integer_offset_x_pix,
        offy=-rectwv_coeff.global_integer_offset_y_pix
    )

    # save output file
    save_ndarray_to_fits(
        array=image2d_flatfielded,
        file_name=args.outfile,
        main_header=header,
        overwrite=True
    )
    print('>>> Saving file ' + args.outfile.name)


if __name__ == "__main__":
    main()
