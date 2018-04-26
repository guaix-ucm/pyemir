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

from numina.array.wavecalib.resample import resample_image2d_flux
from numina.tools.arg_file_is_new import arg_file_is_new

from emirdrp.instrument.dtu_configuration import DtuConfiguration
from emirdrp.products import RectWaveCoeff
from emirdrp.processing.wavecal import Slitlet2D
from emirdrp.processing.wavecal import set_wv_parameters
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers

from .rect_wpoly_for_mos import islitlet_progress
from .save_ndarray_to_fits import save_ndarray_to_fits

from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: apply rectification and wavelength '
                    'calibration polynomials for the CSU configuration of a '
                    'particular image'
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
                        help="Output FITS file with rectified and "
                             "wavelength calibrated image",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))

    # optional arguments
    parser.add_argument("--resampling",
                        help="Resampling method: 1 -> nearest neighbor, "
                             "2 -> linear interpolation (default)",
                        default=2, type=int,
                        choices=(1, 2))
    parser.add_argument("--ignore_DTUconf",
                        help="Ignore DTU configurations differences between "
                             "transformation and input image",
                        action="store_true")
    parser.add_argument("--outfile_rectified_only",
                        help="Output FITS file with rectified image (without "
                             "wavelength calibration!)",
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
    dtu_conf_fitsfile = DtuConfiguration.define_from_fits(args.fitsfile)
    dtu_conf_jsonfile = DtuConfiguration.define_from_dictionary(
        rectwv_coeff.meta_info['dtu_configuration'])
    if dtu_conf_fitsfile != dtu_conf_jsonfile:
        print('DTU configuration (FITS file):\n\t', dtu_conf_fitsfile)
        print('DTU configuration (JSON file):\n\t', dtu_conf_jsonfile)
        if args.ignore_DTUconf:
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

    # ---

    if args.outfile_rectified_only is not None:
        image2d_rectified = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))
        image2d_unrectified = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

        for islitlet in list_valid_islitlets:
            if args.debugplot == 0:
                islitlet_progress(islitlet, EMIR_NBARS)

            # define Slitlet2D object
            slt = Slitlet2D(islitlet=islitlet,
                            rectwv_coeff=rectwv_coeff,
                            debugplot=args.debugplot)

            # minimum and maximum useful scan (pixel in the spatial direction)
            # for the rectified slitlet
            nscan_min, nscan_max = nscan_minmax_frontiers(
                slt.y0_frontier_lower,
                slt.y0_frontier_upper,
                resize=False
            )
            # extract 2D image corresponding to the selected slitlet: note that
            # in this case we are not using select_unrectified_slitlets()
            # because it introduces extra zero pixels in the slitlet frontiers
            slitlet2d = slt.extract_slitlet2d(image2d)

            # rectify image
            slitlet2d_rect = slt.rectify(slitlet2d,
                                         resampling=args.resampling)
            slitlet2d_unrect = slt.rectify(slitlet2d_rect,
                                           resampling=args.resampling,
                                           inverse=True)

            ii1 = nscan_min - slt.bb_ns1_orig
            ii2 = nscan_max - slt.bb_ns1_orig + 1

            j1 = slt.bb_nc1_orig - 1
            j2 = slt.bb_nc2_orig
            i1 = slt.bb_ns1_orig - 1 + ii1
            i2 = i1 + ii2 - ii1

            image2d_rectified[i1:i2, j1:j2] = slitlet2d_rect[ii1:ii2, :]
            image2d_unrectified[i1:i2, j1:j2] = slitlet2d_unrect[ii1:ii2, :]

        save_ndarray_to_fits(
            array=[image2d_rectified, image2d_unrectified],
            file_name=args.outfile_rectified_only,
            cast_to_float=[True] * 2,
            overwrite=True
        )

    # ---

    # relevant wavelength calibration parameters for rectified and wavelength
    # calibrated image
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    crpix1_enlarged = wv_parameters['crpix1_enlarged']
    crval1_enlarged = wv_parameters['crval1_enlarged']
    cdelt1_enlarged = wv_parameters['cdelt1_enlarged']
    naxis1_enlarged = wv_parameters['naxis1_enlarged']

    # initialize rectified and wavelength calibrated image
    image2d_rectified_wv = np.zeros((EMIR_NAXIS2, naxis1_enlarged))

    # main loop
    import time  # TODO: remove this and print(time.ctime()) below
    for islitlet in list_valid_islitlets:
        if args.debugplot == 0:
            if islitlet == list_valid_islitlets[0]:
                print(time.ctime())
            islitlet_progress(islitlet, EMIR_NBARS)
            if islitlet == list_valid_islitlets[-1]:
                print(' ')
                print(time.ctime())

        # define Slitlet2D object
        slt = Slitlet2D(islitlet=islitlet,
                        rectwv_coeff=rectwv_coeff,
                        debugplot=args.debugplot)

        if abs(args.debugplot) >= 10:
            print(slt)

        # extract (distorted) slitlet from the initial image
        slitlet2d = slt.extract_slitlet2d(image2d)

        # rectify slitlet
        slitlet2d_rect = slt.rectify(slitlet2d, resampling=args.resampling)

        # wavelength calibration of the rectifed slitlet
        slitlet2d_rect_wv = resample_image2d_flux(
            image2d_orig=slitlet2d_rect,
            naxis1=naxis1_enlarged,
            cdelt1=cdelt1_enlarged,
            crval1=crval1_enlarged,
            crpix1=crpix1_enlarged,
            coeff=slt.wpoly
        )

        # minimum and maximum useful scan (pixel in the spatial direction)
        # for the rectified slitlet
        nscan_min, nscan_max = nscan_minmax_frontiers(
            slt.y0_frontier_lower,
            slt.y0_frontier_upper,
            resize=False
        )
        ii1 = nscan_min - slt.bb_ns1_orig
        ii2 = nscan_max - slt.bb_ns1_orig + 1
        i1 = slt.bb_ns1_orig - 1 + ii1
        i2 = i1 + ii2 - ii1
        image2d_rectified_wv[i1:i2, :] = slitlet2d_rect_wv[ii1:ii2, :]

        # include scan range in FITS header
        header['sltmin' + str(islitlet).zfill(2)] = i1
        header['sltmax' + str(islitlet).zfill(2)] = i2 - 1

    # modify upper limit of previous slitlet in case of overlapping:
    # note that the overlapped scans have been overwritten with the
    # information from the current slitlet!
    for islitlet in list_valid_islitlets:
        cprevious = 'SLTMAX' + str(islitlet - 1).zfill(2)
        if cprevious in header.keys():
            sltmax_previous = header[cprevious]
            cslitlet = 'SLTMIN' + str(islitlet).zfill(2)
            sltmin_current = header[cslitlet]
            if sltmax_previous >= sltmin_current:
                print('WARNING: ' + cslitlet + '=' +
                      str(sltmin_current).zfill(4) +
                      ' overlaps with ' + cprevious + '=' +
                      str(sltmax_previous).zfill(4) + ' ==> ' + cslitlet +
                      ' set to ' + str(sltmin_current - 1).zfill(4))
                header[cprevious] = sltmin_current - 1

    # update wavelength calibration in FITS header
    header['ctype1'] = 'LAMBDA'
    header['ctype2'] = 'PIXEL'
    header['crpix1'] = (crpix1_enlarged, 'reference pixel')
    header['crpix2'] = (1.0, 'reference pixel')
    header['crval1'] = (crval1_enlarged, 'central wavelength at crpix1')
    header['crval2'] = (1.0, 'central value at crpix2')
    header['cdelt1'] = (cdelt1_enlarged, 'linear dispersion (Angstrom/pixel)')
    header['cdelt2'] = (1.0, 'increment')
    header['cunit1'] = ('angstroms', 'units along axis1')
    header['cunit2'] = ('pixel', 'units along axis2')
    header.remove('cd1_1')
    header.remove('cd1_2')
    header.remove('cd2_1')
    header.remove('cd2_2')
    header.remove('PCD1_1')
    header.remove('PCD1_2')
    header.remove('PCD2_1')
    header.remove('PCD2_2')
    header.remove('PCRPIX1')
    header.remove('PCRPIX2')

    save_ndarray_to_fits(
        array=image2d_rectified_wv,
        file_name=args.outfile,
        main_header=header,
        overwrite=True
    )
    print('>>> Saving file ' + args.outfile.name)


if __name__ == "__main__":
    main()
