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
from numpy.polynomial import Polynomial
import os.path
import sys

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow

from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_NAXIS1


def get_boundaries(bounddict_file, slitlet_number):
    """Read the bounddict json file and return the polynomial boundaries.

    Parameters
    ----------
    bounddict_file : file handler
        File containing the bounddict JSON data.
    slitlet_number : int
        Number of slitlet.

    Returns
    -------
    pol_lower_boundary : numpy polynomial
        Polynomial defining the lower boundary of the slitlet.
    pol_upper_boundary : numpy polynomial
        Polynomial defining the upper boundary of the slitlet.
    xmin_lower : float
        Minimum abscissae for the lower boundary.
    xmax_lower : float
        Maximum abscissae for the lower boundary.
    xmin_upper : float
        Minimum abscissae for the upper boundary.
    xmax_upper : float
        Maximum abscissae for the upper boundary.
    csu_bar_slit_center : float
        CSU bar slit center (in mm)

    """

    bounddict = json.loads(open(bounddict_file.name).read())

    # return values in case the requested slitlet number is not defined
    pol_lower_boundary = None
    pol_upper_boundary = None
    xmin_lower = None
    xmax_lower = None
    xmin_upper = None
    xmax_upper = None
    csu_bar_slit_center = None

    # search the slitlet number in bounddict
    slitlet_label = "slitlet" + str(slitlet_number).zfill(2)
    if slitlet_label in bounddict['contents'].keys():
        list_date_obs = list(bounddict['contents'][slitlet_label].keys())
        list_date_obs.sort()
        num_date_obs = len(list_date_obs)
        if num_date_obs == 1:
            date_obs = list_date_obs[0]
            tmp_dict = bounddict['contents'][slitlet_label][date_obs]
            pol_lower_boundary = Polynomial(tmp_dict['boundary_coef_lower'])
            pol_upper_boundary = Polynomial(tmp_dict['boundary_coef_upper'])
            xmin_lower = tmp_dict['boundary_xmin_lower']
            xmax_lower = tmp_dict['boundary_xmax_lower']
            xmin_upper = tmp_dict['boundary_xmin_upper']
            xmax_upper = tmp_dict['boundary_xmax_upper']
            csu_bar_slit_center = tmp_dict['csu_bar_slit_center']
        else:
            raise ValueError("num_date_obs =", num_date_obs,
                             " (must be 1)")
    else:
        print("WARNING: slitlet number " + str(slitlet_number) +
              " is not available in " + bounddict_file.name)

    # return result
    return pol_lower_boundary, pol_upper_boundary, \
           xmin_lower, xmax_lower, xmin_upper, xmax_upper, \
           csu_bar_slit_center


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("fitsfile",
                        help="FITS file name to be displayed",
                        type=argparse.FileType('rb'))
    parser.add_argument("--bounddict", required=True,
                        help="bounddict file name",
                        type=argparse.FileType('rt'))
    parser.add_argument("--tuple_slit_numbers", required=True,
                        help="Tuple n1[,n2[,step]] to define slitlet numbers")

    # optional arguments
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args()

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # read slitlet numbers to be computed
    tmp_str = args.tuple_slit_numbers.split(",")
    if len(tmp_str) == 3:
        if int(tmp_str[0]) < 1:
            raise ValueError("Invalid slitlet number < 1")
        if int(tmp_str[1]) > EMIR_NBARS:
            raise ValueError("Invalid slitlet number > EMIR_NBARS")
        list_slitlets = range(int(tmp_str[0]),
                              int(tmp_str[1])+1,
                              int(tmp_str[2]))
    elif len(tmp_str) == 2:
        if int(tmp_str[0]) < 1:
            raise ValueError("Invalid slitlet number < 1")
        if int(tmp_str[1]) > EMIR_NBARS:
            raise ValueError("Invalid slitlet number > EMIR_NBARS")
        list_slitlets = range(int(tmp_str[0]),
                              int(tmp_str[1])+1,
                              1)
    elif len(tmp_str) == 1:
        if int(tmp_str[0]) < 1:
            raise ValueError("Invalid slitlet number < 1")
        if int(tmp_str[0]) > EMIR_NBARS:
            raise ValueError("Invalid slitlet number > EMIR_NBARS")
        list_slitlets = [int(tmp_str[0])]
    else:
        raise ValueError("Invalid tuple for slitlet numbers")

    # read input FITS file
    hdulist = fits.open(args.fitsfile.name)
    image_header = hdulist[0].header
    image2d = hdulist[0].data
    hdulist.close()

    naxis1 = image_header['naxis1']
    naxis2 = image_header['naxis2']

    if image2d.shape != (naxis2, naxis1):
        raise ValueError("Unexpected error with NAXIS1, NAXIS2")

    # remove path from fitsfile
    sfitsfile = os.path.basename(args.fitsfile.name)

    # check that the FITS file has been obtained with EMIR
    instrument = image_header['instrume']
    if instrument != 'EMIR':
        raise ValueError("INSTRUME keyword is not 'EMIR'!")

    # read GRISM, FILTER and ROTANG from FITS header
    grism = image_header['grism']
    spfilter = image_header['filter']
    rotang = image_header['rotang']

    # display full image
    ax = ximshow(image2d=image2d,
                 title=sfitsfile + "\ngrism=" + grism +
                       ", filter=" + spfilter +
                       ", rotang=" + str(round(rotang, 2)),
                 image_bbox=(1, naxis1, 1, naxis2), show=False)

    # overplot boundaries for each slitlet
    for slitlet_number in list_slitlets:
        pol_lower_boundary, pol_upper_boundary, \
        xmin_lower, xmax_lower, xmin_upper, xmax_upper,  \
        csu_bar_slit_center = \
            get_boundaries(args.bounddict, slitlet_number)
        if (pol_lower_boundary is not None) and \
                (pol_upper_boundary is not None):
            xp = np.linspace(start=xmin_lower, stop=xmax_lower, num=1000)
            yp = pol_lower_boundary(xp)
            ax.plot(xp, yp, 'g-')
            xp = np.linspace(start=xmin_upper, stop=xmax_upper, num=1000)
            yp = pol_upper_boundary(xp)
            ax.plot(xp, yp, 'b-')
            # slitlet label
            yc_lower = pol_lower_boundary(EMIR_NAXIS1 / 2 + 0.5)
            yc_upper = pol_upper_boundary(EMIR_NAXIS1 / 2 + 0.5)
            tmpcolor = ['r', 'b'][slitlet_number % 2]
            xcsu = EMIR_NAXIS1 * csu_bar_slit_center/341.5
            ax.text(xcsu, (yc_lower + yc_upper) / 2,
                    str(slitlet_number),
                    fontsize=10, va='center', ha='center',
                    bbox=dict(boxstyle="round,pad=0.1",
                              fc="white", ec="grey"),
                    color=tmpcolor, fontweight='bold',
                    backgroundcolor='white')

    # show plot
    pause_debugplot(12, pltshow=True)


if __name__ == "__main__":

    main()
