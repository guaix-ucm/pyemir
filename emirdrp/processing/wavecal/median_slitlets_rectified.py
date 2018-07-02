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

from emirdrp.instrument.csu_configuration import CsuConfiguration

from numina.array.display.ximshow import ximshow
from numina.array.wavecalib.fix_pix_borders import define_mask_borders
from numina.tools.arg_file_is_new import arg_file_is_new

from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_NPIXPERSLIT_RECTIFIED

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_MINIMUM_SLITLET_WIDTH_MM
from emirdrp.core import EMIR_MAXIMUM_SLITLET_WIDTH_MM


def median_slitlets_rectified(
        input_image,
        mode=0,
        minimum_slitlet_width_mm=EMIR_MINIMUM_SLITLET_WIDTH_MM,
        maximum_slitlet_width_mm=EMIR_MAXIMUM_SLITLET_WIDTH_MM,
        debugplot=0
    ):
    """Compute median spectrum for each slitlet.

    Parameters
    ----------
    input_image : HDUList object
        Input 2D image.
    mode : int
        Indicate desired result:
        0 : image with the same size as the input image, with the
            median spectrum of each slitlet spanning all the spectra
            of the corresponding slitlet
        1 : image with 55 spectra, containing the median spectra of
            each slitlet
        2 : single collapsed median spectrum, using exclusively the
            useful pixels from the input image
    minimum_slitlet_width_mm : float
        Minimum slitlet width (mm) for a valid slitlet.
    maximum_slitlet_width_mm : float
        Maximum slitlet width (mm) for a valid slitlet.
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed. The valid codes are defined in
        numina.array.display.pause_debugplot.

    Returns
    -------
    image_median : HDUList object
        Output image.

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
    if mode == 0:
        image2d_median = np.zeros((naxis2, naxis1))
    else:
        image2d_median = np.zeros((EMIR_NBARS, naxis1))

    # main loop
    for i in range(EMIR_NBARS):
        ns1 = i * EMIR_NPIXPERSLIT_RECTIFIED + 1
        ns2 = ns1 + EMIR_NPIXPERSLIT_RECTIFIED - 1
        sp_median = np.median(image2d[(ns1-1):ns2, :], axis=0)

        if mode == 0:
            image2d_median[(ns1-1):ns2, :] = np.tile(
                sp_median, (EMIR_NPIXPERSLIT_RECTIFIED, 1)
            )
        else:
            image2d_median[i] = np.copy(sp_median)

    if mode == 2:
        # get CSU configuration from FITS header
        csu_config = CsuConfiguration.define_from_header(image_header)

        # define wavelength calibration parameters
        crpix1 = image_header['crpix1']
        crval1 = image_header['crval1']
        cdelt1 = image_header['cdelt1']

        # segregate slitlets
        list_useful_slitlets = csu_config.widths_in_range_mm(
            minwidth=minimum_slitlet_width_mm,
            maxwidth=maximum_slitlet_width_mm
        )
        list_not_useful_slitlets = [i for i in list(range(1, EMIR_NBARS + 1))
                                    if i not in list_useful_slitlets]
        if abs(debugplot) != 0:
            print('>>> list_useful_slitlets....:', list_useful_slitlets)
            print('>>> list_not_useful_slitlets:', list_not_useful_slitlets)

        # define mask from array data
        mask2d, borders = define_mask_borders(image2d_median, sought_value=0)
        if abs(debugplot) % 10 != 0:
            ximshow(mask2d.astype(int), z1z2=(-.2, 1.2), crpix1=crpix1,
                    crval1=crval1, cdelt1=cdelt1, debugplot=debugplot)

        # update mask with unused slitlets
        for islitlet in list_not_useful_slitlets:
            mask2d[islitlet - 1, :] = np.array([True] * naxis1)
        if abs(debugplot) % 10 != 0:
            ximshow(mask2d.astype(int), z1z2=(-.2, 1.2), crpix1=crpix1,
                    crval1=crval1, cdelt1=cdelt1, debugplot=debugplot)

        # useful image pixels
        image2d_masked = image2d_median * (1 - mask2d.astype(int))
        if abs(debugplot) % 10 != 0:
            ximshow(image2d_masked, crpix1=crpix1, crval1=crval1,
                    cdelt1=cdelt1, debugplot=debugplot)

        # masked image
        image2d_masked = np.ma.masked_array(image2d_median, mask=mask2d)
        # median spectrum
        image1d_median = np.ma.median(image2d_masked, axis=0).data

        image_median = fits.PrimaryHDU(data=image1d_median,
                                       header=image_header)

    else:
        image_median = fits.PrimaryHDU(data=image2d_median,
                                       header=image_header)

    return fits.HDUList([image_median])


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: compute median spectrum for each slitlet'
    )

    # positional arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file name",
                        type=argparse.FileType('rb'))
    parser.add_argument("outfile",
                        help="Output FITS file name",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))

    # optional arguments
    parser.add_argument("--mode",
                        help="Output type: 0 -> full frame (default), "
                             "1 -> individual slitlets, "
                             "2 -> collapsed single spectrum)",
                        default=0, type=int,
                        choices=[0, 1, 2])
    parser.add_argument("--minimum_slitlet_width_mm",
                        help="Minimum slitlet width (mm) for --mode 2 "
                             "(default=0)",
                        default=EMIR_MINIMUM_SLITLET_WIDTH_MM, type=float)
    parser.add_argument("--maximum_slitlet_width_mm",
                        help="Maximum slitlet width (mm) for --mode 2 "
                             "(default=" +
                             str(EMIR_MAXIMUM_SLITLET_WIDTH_MM) + ")",
                        default=EMIR_MAXIMUM_SLITLET_WIDTH_MM, type=float)
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

    # read input FITS file
    hdulist = fits.open(args.fitsfile)

    image_median = median_slitlets_rectified(
        hdulist,
        mode=args.mode,
        minimum_slitlet_width_mm=args.minimum_slitlet_width_mm,
        maximum_slitlet_width_mm=args.maximum_slitlet_width_mm,
        debugplot=args.debugplot
    )

    # save result
    image_median.writeto(args.outfile, overwrite=True)


if __name__ == "__main__":

    main()
