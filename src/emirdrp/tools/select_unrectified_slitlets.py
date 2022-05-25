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
import json
import numpy as np
import os.path
import sys

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow
from numina.tools.arg_file_is_new import arg_file_is_new
from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers

from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.fit_boundaries import expected_distorted_frontiers
from emirdrp.tools.fit_boundaries import overplot_boundaries_from_params
from emirdrp.tools.fit_boundaries import overplot_frontiers_from_params
from emirdrp.tools.list_slitlets_from_string import list_slitlets_from_string

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2


def select_unrectified_slitlet(image2d, islitlet, csu_bar_slit_center,
                               params, parmodel, maskonly):
    """Returns image with the indicated slitlet (zero anywhere else).

    Parameters
    ----------
    image2d : numpy array
        Initial image from which the slitlet data will be extracted.
    islitlet : int
        Slitlet number.
    csu_bar_slit_center : float
        CSU bar slit center.
    params : :class:`~lmfit.parameter.Parameters`
        Parameters to be employed in the prediction of the distorted
        boundaries.
    parmodel : str
        Model to be assumed. Allowed values are 'longslit' and
        'multislit'.
    maskonly : bool
        If True, returns simply a mask (1 in the slitlet region and
        zero anywhere else.

    Returns
    -------
    image2d_output : numpy array
        2D image with the pixel information corresponding to the
        selected slitlet and zero everywhere else.

    """

    # protection
    if image2d.shape != (EMIR_NAXIS2, EMIR_NAXIS1):
        raise ValueError("NAXIS1, NAXIS2 unexpected for EMIR detector")

    # initialize image output
    image2d_output = np.zeros_like(image2d)

    # expected slitlet frontiers
    list_expected_frontiers = expected_distorted_frontiers(
        islitlet, csu_bar_slit_center,
        params, parmodel, numpts=101, deg=5, debugplot=0
    )
    pol_lower_expected = list_expected_frontiers[0].poly_funct
    pol_upper_expected = list_expected_frontiers[1].poly_funct

    # main loop: compute for each channel the minimum and maximum scan
    for j in range(EMIR_NAXIS1):
        xchannel = j + 1
        y0_lower = pol_lower_expected(xchannel)
        y0_upper = pol_upper_expected(xchannel)
        n1, n2 = nscan_minmax_frontiers(y0_frontier_lower=y0_lower,
                                        y0_frontier_upper=y0_upper,
                                        resize=True)
        # note that n1 and n2 are scans (ranging from 1 to NAXIS2)
        if maskonly:
            image2d_output[(n1 - 1):n2, j] = np.repeat(
                [1.0], (n2 - n1 + 1)
            )
        else:
            image2d_output[(n1 - 1):n2, j] = image2d[(n1 - 1):n2, j]

    return image2d_output


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("fitsfile",
                        help="FITS file name to be displayed",
                        type=argparse.FileType('rb'))
    parser.add_argument("--fitted_bound_param", required=True,
                        help="JSON file with fitted boundary coefficients "
                             "corresponding to the multislit model",
                        type=argparse.FileType('rt'))
    parser.add_argument("--slitlets", required=True,
                        help="Slitlet selection: string between double "
                             "quotes providing tuples of the form "
                             "n1[,n2[,step]]",
                        type=str)

    # optional arguments
    parser.add_argument("--outfile",
                        help="Output FITS file name",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))
    parser.add_argument("--maskonly",
                        help="Generate mask for the indicated slitlets",
                        action="store_true")
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting/debugging" +
                             " (default=0)",
                        type=int, default=0,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args()

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # read input FITS file
    hdulist_image = fits.open(args.fitsfile.name)
    image_header = hdulist_image[0].header
    image2d = hdulist_image[0].data

    naxis1 = image_header['naxis1']
    naxis2 = image_header['naxis2']

    if image2d.shape != (naxis2, naxis1):
        raise ValueError("Unexpected error with NAXIS1, NAXIS2")

    if image2d.shape != (EMIR_NAXIS2, EMIR_NAXIS1):
        raise ValueError("NAXIS1, NAXIS2 unexpected for EMIR detector")

    # remove path from fitsfile
    if args.outfile is None:
        sfitsfile = os.path.basename(args.fitsfile.name)
    else:
        sfitsfile = os.path.basename(args.outfile.name)

    # check that the FITS file has been obtained with EMIR
    instrument = image_header['instrume']
    if instrument != 'EMIR':
        raise ValueError("INSTRUME keyword is not 'EMIR'!")

    # read GRISM, FILTER and ROTANG from FITS header
    grism = image_header['grism']
    spfilter = image_header['filter']
    rotang = image_header['rotang']

    # read fitted_bound_param JSON file
    fittedpar_dict = json.loads(open(args.fitted_bound_param.name).read())
    params = bound_params_from_dict(fittedpar_dict)
    if abs(args.debugplot) in [21, 22]:
        params.pretty_print()

    parmodel = fittedpar_dict['meta_info']['parmodel']
    if parmodel != 'multislit':
        raise ValueError("Unexpected parameter model: ", parmodel)

    # define slitlet range
    islitlet_min = fittedpar_dict['tags']['islitlet_min']
    islitlet_max = fittedpar_dict['tags']['islitlet_max']
    list_islitlet = list_slitlets_from_string(
        s=args.slitlets,
        islitlet_min=islitlet_min,
        islitlet_max=islitlet_max
    )

    # read CsuConfiguration object from FITS file
    csu_config = CsuConfiguration.define_from_fits(args.fitsfile)

    # define csu_bar_slit_center associated to each slitlet
    list_csu_bar_slit_center = []
    for islitlet in list_islitlet:
        list_csu_bar_slit_center.append(
            csu_config.csu_bar_slit_center(islitlet))

    # initialize output data array
    image2d_output = np.zeros((naxis2, naxis1))

    # main loop
    for islitlet, csu_bar_slit_center in \
            zip(list_islitlet, list_csu_bar_slit_center):
        image2d_tmp = select_unrectified_slitlet(
            image2d=image2d,
            islitlet=islitlet,
            csu_bar_slit_center=csu_bar_slit_center,
            params=params,
            parmodel=parmodel,
            maskonly=args.maskonly
        )
        image2d_output += image2d_tmp

    # update the array of the output file
    hdulist_image[0].data = image2d_output

    # save output FITS file
    hdulist_image.writeto(args.outfile)

    # close original image
    hdulist_image.close()

    # display full image
    if abs(args.debugplot) % 10 != 0:
        ax = ximshow(image2d=image2d_output,
                     title=sfitsfile + "\n" + args.slitlets,
                     image_bbox=(1, naxis1, 1, naxis2), show=False)

        # overplot boundaries
        overplot_boundaries_from_params(
            ax=ax,
            params=params,
            parmodel=parmodel,
            list_islitlet=list_islitlet,
            list_csu_bar_slit_center=list_csu_bar_slit_center
        )

        # overplot frontiers
        overplot_frontiers_from_params(
            ax=ax,
            params=params,
            parmodel=parmodel,
            list_islitlet=list_islitlet,
            list_csu_bar_slit_center=list_csu_bar_slit_center,
            micolors=('b', 'b'), linetype='-',
            labels=False    # already displayed with the boundaries
        )

        # show plot
        pause_debugplot(12, pltshow=True)


if __name__ == "__main__":

    main()
