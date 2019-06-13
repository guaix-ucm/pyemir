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
from sklearn.cluster import KMeans
import sys

from numina.array.display.matplotlib_qt import plt
from numina.array.display.matplotlib_qt import patches as patches
from numina.array.display.matplotlib_qt import set_window_geometry
from numina.array.display.pause_debugplot import pause_debugplot
from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.processing.wavecal.set_wv_parameters import set_wv_parameters

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NBARS
from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def display_slitlet_arrangement(fileobj,
                                grism=None,
                                spfilter=None,
                                longslits=False,
                                bbox=None,
                                adjust=None,
                                geometry=None,
                                debugplot=0):
    """Display slitlet arrangment from CSUP keywords in FITS header.

    Parameters
    ----------
    fileobj : file object
        FITS or TXT file object.
    grism : str
        Grism name.
    spfilter : str
        Filter name.
    longslits : bool
        If True, display (pseudo) longslits.
    bbox : tuple of 4 floats
        If not None, values for xmin, xmax, ymin and ymax.
    adjust : bool
        Adjust X range according to minimum and maximum csu_bar_left
        and csu_bar_right (note that this option is overriden by 'bbox')
    geometry : tuple (4 integers) or None
        x, y, dx, dy values employed to set the Qt backend geometry.
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed. The valid codes are defined in
        numina.array.display.pause_debugplot

    Returns
    -------
    csu_bar_left : list of floats
        Location (mm) of the left bar for each slitlet.
    csu_bar_right : list of floats
        Location (mm) of the right bar for each slitlet, using the
        same origin employed for csu_bar_left (which is not the
        value stored in the FITS keywords.
    csu_bar_slit_center : list of floats
        Middle point (mm) in between the two bars defining a slitlet.
    csu_bar_slit_width : list of floats
        Slitlet width (mm), computed as the distance between the two
        bars defining the slitlet.

    """

    if fileobj.name[-4:] == ".txt":
        if grism is None:
            raise ValueError("Undefined grism!")
        if spfilter is None:
            raise ValueError("Undefined filter!")
        # define CsuConfiguration object
        csu_config = CsuConfiguration()
        csu_config._csu_bar_left = []
        csu_config._csu_bar_right = []
        csu_config._csu_bar_slit_center = []
        csu_config._csu_bar_slit_width = []

        # since the input filename has been opened with argparse in binary
        # mode, it is necessary to close it and open it in text mode
        fileobj.close()
        # read TXT file
        with open(fileobj.name, mode='rt') as f:
            file_content = f.read().splitlines()
        next_id_bar = 1
        for line in file_content:
            if len(line) > 0:
                if line[0] not in ['#']:
                    line_contents = line.split()
                    id_bar = int(line_contents[0])
                    position = float(line_contents[1])
                    if id_bar == next_id_bar:
                        if id_bar <= EMIR_NBARS:
                            csu_config._csu_bar_left.append(position)
                            next_id_bar = id_bar + EMIR_NBARS
                        else:
                            csu_config._csu_bar_right.append(341.5 - position)
                            next_id_bar = id_bar - EMIR_NBARS + 1
                    else:
                        raise ValueError("Unexpected id_bar:" + str(id_bar))

        # compute slit width and center
        for i in range(EMIR_NBARS):
            csu_config._csu_bar_slit_center.append(
                (csu_config._csu_bar_left[i] + csu_config._csu_bar_right[i])/2
            )
            csu_config._csu_bar_slit_width.append(
                csu_config._csu_bar_right[i] - csu_config._csu_bar_left[i]
            )

    else:
        # read input FITS file
        hdulist = fits.open(fileobj.name)
        image_header = hdulist[0].header
        hdulist.close()

        # additional info from header
        grism = image_header['grism']
        spfilter = image_header['filter']

        # define slitlet arrangement
        csu_config = CsuConfiguration.define_from_fits(fileobj)
        if longslits:
            csu_config.display_pseudo_longslits()
            pause_debugplot(debugplot)

    # determine calibration
    if grism in ["J", "OPEN"] and spfilter == "J":
        wv_parameters = set_wv_parameters("J", "J")
    elif grism in ["H", "OPEN"] and spfilter == "H":
        wv_parameters = set_wv_parameters("H", "H")
    elif grism in ["K", "OPEN"] and spfilter == "Ksp":
        wv_parameters = set_wv_parameters("Ksp", "K")
    elif grism in ["LR", "OPEN"] and spfilter == "YJ":
        wv_parameters = set_wv_parameters("YJ", "LR")
    elif grism in ["LR", "OPEN"] and spfilter == "HK":
        wv_parameters = set_wv_parameters("HK", "LR")
    else:
        raise ValueError("Invalid grism + filter configuration")

    crval1 = wv_parameters['poly_crval1_linear']
    cdelt1 = wv_parameters['poly_cdelt1_linear']

    wvmin_useful = wv_parameters['wvmin_useful']
    wvmax_useful = wv_parameters['wvmax_useful']

    # display arrangement
    if abs(debugplot) >= 10:
        print("slit     left    right   center   width   min.wave   max.wave")
        print("====  =======  =======  =======   =====   ========   ========")
        for i in range(EMIR_NBARS):
            ibar = i + 1
            csu_crval1 = crval1(csu_config.csu_bar_slit_center(ibar))
            csu_cdelt1 = cdelt1(csu_config.csu_bar_slit_center(ibar))
            csu_crvaln = csu_crval1 + (EMIR_NAXIS1 - 1) * csu_cdelt1
            if wvmin_useful is not None:
                csu_crval1 = np.amax([csu_crval1, wvmin_useful])
            if wvmax_useful is not None:
                csu_crvaln = np.amin([csu_crvaln, wvmax_useful])
            print("{0:4d} {1:8.3f} {2:8.3f} {3:8.3f} {4:7.3f}   "
                  "{5:8.2f}   {6:8.2f}".format(
                ibar, csu_config.csu_bar_left(ibar),
                csu_config.csu_bar_right(ibar),
                csu_config.csu_bar_slit_center(ibar),
                csu_config.csu_bar_slit_width(ibar),
                csu_crval1, csu_crvaln)
            )
        print(
            "---> {0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} <- mean (all)".format(
                np.mean(csu_config._csu_bar_left),
                np.mean(csu_config._csu_bar_right),
                np.mean(csu_config._csu_bar_slit_center),
                np.mean(csu_config._csu_bar_slit_width)
            )
        )
        print(
            "---> {0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} <- mean (odd)".format(
                np.mean(csu_config._csu_bar_left[::2]),
                np.mean(csu_config._csu_bar_right[::2]),
                np.mean(csu_config._csu_bar_slit_center[::2]),
                np.mean(csu_config._csu_bar_slit_width[::2])
            )
        )
        print(
            "---> {0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} <- mean (even)".format(
                np.mean(csu_config._csu_bar_left[1::2]),
                np.mean(csu_config._csu_bar_right[1::2]),
                np.mean(csu_config._csu_bar_slit_center[1::2]),
                np.mean(csu_config._csu_bar_slit_width[1::2])
            )
        )

    # display slit arrangement
    if abs(debugplot) % 10 != 0:
        fig = plt.figure()
        set_window_geometry(geometry)
        ax = fig.add_subplot(111)
        if bbox is None:
            if adjust:
                xmin = min(csu_config._csu_bar_left)
                xmax = max(csu_config._csu_bar_right)
                dx = xmax - xmin
                if dx == 0:
                    dx = 1
                xmin -= dx/20
                xmax += dx/20
                ax.set_xlim(xmin, xmax)
            else:
                ax.set_xlim(0., 341.5)
            ax.set_ylim(0, 56)
        else:
            ax.set_xlim(bbox[0], bbox[1])
            ax.set_ylim(bbox[2], bbox[3])
        ax.set_xlabel('csu_bar_position (mm)')
        ax.set_ylabel('slit number')
        for i in range(EMIR_NBARS):
            ibar = i + 1
            ax.add_patch(patches.Rectangle(
                (csu_config.csu_bar_left(ibar), ibar-0.5),
                csu_config.csu_bar_slit_width(ibar), 1.0))
            ax.plot([0., csu_config.csu_bar_left(ibar)], [ibar, ibar],
                    '-', color='gray')
            ax.plot([csu_config.csu_bar_right(ibar), 341.5],
                    [ibar, ibar], '-', color='gray')
        plt.title("File: " + fileobj.name + "\ngrism=" + grism +
                  ", filter=" + spfilter)
        pause_debugplot(debugplot, pltshow=True)

    # return results
    return csu_config._csu_bar_left, csu_config._csu_bar_right, \
           csu_config._csu_bar_slit_center, csu_config._csu_bar_slit_width


def display_slitlet_histogram(csu_bar_slit_width,
                              n_clusters=2,
                              geometry=None,
                              debugplot=0):
    """

    Find separations between groups of slitlet widths.

    Parameters
    ----------
    csu_bar_slit_width : numpy array
        Array containing the csu_bar_slit_center values.
    n_clusters : int
        Number of slitlet groups to be sought.
    geometry : tuple (4 integers) or None
        x, y, dx, dy values employed to set the Qt backend geometry.
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed. The valid codes are defined in
        numina.array.display.pause_debugplot

    Returns
    -------
    TBD

    """

    # protections
    if n_clusters < 2:
        raise ValueError("n_clusters must be >= 2")

    # determine useful slitlets from their widths
    km = KMeans(n_clusters=n_clusters)
    fit = km.fit(np.array(csu_bar_slit_width).reshape(-1, 1))
    separator_list = []
    for i in range(n_clusters - 1):
        peak1 = fit.cluster_centers_[i][0]
        peak2 = fit.cluster_centers_[i+1][0]
        separator = (peak1 + peak2) / 2
        separator_list.append(separator)
        array_csu_bar_slit_width = np.array(csu_bar_slit_width)
        list_ok = []
        list_not_ok = []
        for k in range(EMIR_NBARS):
            if array_csu_bar_slit_width[k] > separator:
                list_ok.append(k + 1)
            else:
                list_not_ok.append(k + 1)
        print('\nNumber of slitlets with width > separator: {}'.format(
            sum(np.array(csu_bar_slit_width) > separator)
        ))
        print(list_ok)
        print('\nNumber of slitlets with width < separator: {}'.format(
            sum(np.array(csu_bar_slit_width) < separator)
        ))
        print(list_not_ok)
        print('\n--->  separator: {0:7.3f}'.format(separator))

    # display histogram
    if abs(debugplot) % 10 != 0:
        fig = plt.figure()
        set_window_geometry(geometry)
        ax = fig.add_subplot(111)
        ax.hist(csu_bar_slit_width, bins=100)
        ax.set_xlabel('slit width (mm)')
        ax.set_ylabel('number of slitlets')
        for separator in separator_list:
            ax.axvline(separator, color='C1', linestyle='dashed')
        pause_debugplot(debugplot, pltshow=True)


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: display arrangement of EMIR CSU bars'
    )

    # positional arguments
    parser.add_argument("filename",
                        help="FITS files (wildcards accepted) or single TXT "
                             "file with CSU configuration from OSP",
                        type=argparse.FileType('rb'),
                        nargs='+')

    # optional arguments
    parser.add_argument("--grism",
                        help="Grism (J, H, K, LR)",
                        choices=["J", "H", "K", "LR"])
    parser.add_argument("--filter",
                        help="Filter (J, H, Ksp, YJ, HK)",
                        choices=["J", "H", "Ksp", "YJ", "HK"])
    parser.add_argument("--n_clusters",
                        help="Display histogram of slitlet widths",
                        default=0, type=int)
    parser.add_argument("--longslits",
                        help="Display (pseudo) longslits built with "
                             "consecutive slitlets",
                        action="store_true")
    parser.add_argument("--bbox",
                        help="Bounding box tuple xmin,xmax,ymin,ymax "
                             "indicating plot limits")
    parser.add_argument("--adjust",
                        help="Adjust X range according to minimum and maximum"
                             " csu_bar_left and csu_bar_right (note that this "
                             "option is overriden by --bbox",
                        action='store_true')
    parser.add_argument("--geometry",
                        help="Tuple x,y,dx,dy indicating window geometry",
                        default="0,0,640,480")
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting & debugging options"
                             " (default=12)",
                        default=12, type=int,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")
    args = parser.parse_args(args)

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

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

    # read bounding box
    if args.bbox is None:
        bbox = None
    else:
        str_bbox = args.bbox.split(",")
        xmin, xmax, ymin, ymax = [int(str_bbox[i]) for i in range(4)]
        bbox = xmin, xmax, ymin, ymax

    # if input file is a single txt file, assume it is a list of FITS files
    if len(args.filename) == 1:
        list_fits_file_objects = [args.filename[0]]
    else:
        list_fits_file_objects = args.filename

    # total number of files to be examined
    nfiles = len(list_fits_file_objects)

    # declare arrays to store CSU values
    csu_bar_left = np.zeros((nfiles, EMIR_NBARS))
    csu_bar_right = np.zeros((nfiles, EMIR_NBARS))
    csu_bar_slit_center = np.zeros((nfiles, EMIR_NBARS))
    csu_bar_slit_width = np.zeros((nfiles, EMIR_NBARS))

    # display CSU bar arrangement
    for ifile, fileobj in enumerate(list_fits_file_objects):
        print("\nFile " + str(ifile+1) + "/" + str(nfiles) + ": " +
              fileobj.name)
        csu_bar_left[ifile, :], csu_bar_right[ifile, :], \
        csu_bar_slit_center[ifile, :], csu_bar_slit_width[ifile, :] = \
            display_slitlet_arrangement(
                fileobj,
                grism=args.grism,
                spfilter=args.filter,
                longslits=args.longslits,
                bbox=bbox,
                adjust=args.adjust,
                geometry=geometry,
                debugplot=args.debugplot
            )
        if args.n_clusters >= 2:
            display_slitlet_histogram(
                csu_bar_slit_width[ifile, :],
                n_clusters=args.n_clusters,
                geometry=geometry,
                debugplot=args.debugplot
            )

    # print summary of comparison between files
    if nfiles > 1:
        std_csu_bar_left = np.zeros(EMIR_NBARS)
        std_csu_bar_right = np.zeros(EMIR_NBARS)
        std_csu_bar_slit_center = np.zeros(EMIR_NBARS)
        std_csu_bar_slit_width = np.zeros(EMIR_NBARS)
        if args.debugplot >= 10:
            print("\n   STANDARD DEVIATION BETWEEN IMAGES")
            print("slit     left    right   center   width")
            print("====  =======  =======  =======   =====")
            for i in range(EMIR_NBARS):
                ibar = i + 1
                std_csu_bar_left[i] = np.std(csu_bar_left[:, i])
                std_csu_bar_right[i] = np.std(csu_bar_right[:, i])
                std_csu_bar_slit_center[i] = np.std(csu_bar_slit_center[:, i])
                std_csu_bar_slit_width[i] = np.std(csu_bar_slit_width[:, i])
                print("{0:4d} {1:8.3f} {2:8.3f} {3:8.3f} {4:7.3f}".format(
                    ibar,
                    std_csu_bar_left[i],
                    std_csu_bar_right[i],
                    std_csu_bar_slit_center[i],
                    std_csu_bar_slit_width[i]))
            print("====  =======  =======  =======   =====")
            print("MIN: {0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f}".format(
                std_csu_bar_left.min(),
                std_csu_bar_right.min(),
                std_csu_bar_slit_center.min(),
                std_csu_bar_slit_width.min()))
            print("MAX: {0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f}".format(
                std_csu_bar_left.max(),
                std_csu_bar_right.max(),
                std_csu_bar_slit_center.max(),
                std_csu_bar_slit_width.max()))
            print("====  =======  =======  =======   =====")
            print("Total number of files examined:", nfiles)

    # stop program execution
    if len(list_fits_file_objects) > 1:
        pause_debugplot(12, optional_prompt="Press RETURN to STOP")


if __name__ == "__main__":
    main()
