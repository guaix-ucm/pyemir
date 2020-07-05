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

from emirdrp.instrument.dtuconf import DtuConf
from emirdrp.products import RectWaveCoeff
from emirdrp.processing.wavecal.slitlet2d import Slitlet2D
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers

from .rect_wpoly_for_mos import islitlet_progress
from .save_ndarray_to_fits import save_ndarray_to_fits

from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NPIXPERSLIT_RECTIFIED

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: apply rectification polynomials '
                    'for the CSU configuration of a particular image'
    )

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file",
                        type=argparse.FileType('rb'))
    parser.add_argument("--rectwv_coeff", required=True,
                        help="Input JSON file with rectification and "
                             "wavelength calibration coefficients",
                        type=argparse.FileType('rt'))
    parser.add_argument("--outfile", required=True,
                        help="Output FITS file with rectified image",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))

    # optional arguments
    parser.add_argument("--resampling",
                        help="Resampling method: 1 -> nearest neighbor, "
                             "2 -> linear interpolation (default)",
                        default=2, type=int,
                        choices=(1, 2))
    parser.add_argument("--ignore_dtu_configuration",
                        help="Ignore DTU configurations differences between "
                             "transformation and input image",
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
    rectwv_coeff = RectWaveCoeff._datatype_load(
        args.rectwv_coeff.name)

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
    if filter_name != rectwv_coeff.tags['filter']:
        raise ValueError("Filter name does not match!")
    grism_name = header['grism']
    if grism_name != rectwv_coeff.tags['grism']:
        raise ValueError("Filter name does not match!")
    if abs(args.debugplot) >= 10:
        print('>>> grism.......:', grism_name)
        print('>>> filter......:', filter_name)

    # check that the DTU configurations are compatible
    with fits.open(args.fitsfile) as hdulist:
        dtu_conf_fitsfile = DtuConf.from_img(hdulist)
    dtu_conf_jsonfile = DtuConf.from_values(
        **rectwv_coeff.meta_info['dtu_configuration']
    )
    if dtu_conf_fitsfile != dtu_conf_jsonfile:
        print('DTU configuration (FITS file):\n\t', dtu_conf_fitsfile)
        print('DTU configuration (JSON file):\n\t', dtu_conf_jsonfile)
        if args.ignore_dtu_configuration:
            print('WARNING: DTU configuration differences found!')
        else:
            raise ValueError("DTU configurations do not match!")
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

    naxis2_enlarged = EMIR_NBARS * EMIR_NPIXPERSLIT_RECTIFIED
    image2d_rectified = np.zeros((naxis2_enlarged, EMIR_NAXIS1))
    image2d_unrectified = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

    for islitlet in list(range(1, EMIR_NBARS + 1)):
        if islitlet in list_valid_islitlets:
            if args.debugplot == 0:
                islitlet_progress(islitlet, EMIR_NBARS, ignore=False)

            # define Slitlet2D object
            slt = Slitlet2D(islitlet=islitlet,
                            rectwv_coeff=rectwv_coeff,
                            debugplot=args.debugplot)

            # extract 2D image corresponding to the selected slitlet: note that
            # in this case we are not using select_unrectified_slitlets()
            # because it introduces extra zero pixels in the slitlet frontiers
            slitlet2d = slt.extract_slitlet2d(image2d)

            # rectify image
            slitlet2d_rect = slt.rectify(slitlet2d,
                                         resampling=args.resampling)

            # minimum and maximum useful row in the full 2d rectified image
            # (starting from 0)
            i1 = slt.iminslt - 1
            i2 = slt.imaxslt

            # minimum and maximum scan in the rectified slitlet
            # (in pixels, from 1 to NAXIS2)
            ii1 = slt.min_row_rectified
            ii2 = slt.max_row_rectified + 1

            # save rectified slitlet in its corresponding location within
            # the full 2d rectified image
            image2d_rectified[i1:i2, :] = slitlet2d_rect[ii1:ii2, :]

            # ---

            # unrectify image
            slitlet2d_unrect = slt.rectify(slitlet2d_rect,
                                           resampling=args.resampling,
                                           inverse=True)

            # minimum and maximum useful scan (pixel in the spatial direction)
            # for the rectified slitlet
            nscan_min, nscan_max = nscan_minmax_frontiers(
                slt.y0_frontier_lower,
                slt.y0_frontier_upper,
                resize=False
            )
            ii1 = nscan_min - slt.bb_ns1_orig
            ii2 = nscan_max - slt.bb_ns1_orig + 1

            j1 = slt.bb_nc1_orig - 1
            j2 = slt.bb_nc2_orig
            i1 = slt.bb_ns1_orig - 1 + ii1
            i2 = i1 + ii2 - ii1

            image2d_unrectified[i1:i2, j1:j2] = slitlet2d_unrect[ii1:ii2, :]

        else:
            if args.debugplot == 0:
                islitlet_progress(islitlet, EMIR_NBARS, ignore=True)

    if args.debugplot == 0:
        print('OK!')

    save_ndarray_to_fits(
        array=[image2d_rectified, image2d_unrectified],
        file_name=args.outfile,
        cast_to_float=[True] * 2,
        overwrite=True
    )
    print('>>> Saving file ' + args.outfile.name)


if __name__ == "__main__":
    main()
