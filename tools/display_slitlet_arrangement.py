from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
import numpy as np
import os.path
import sys

from numina.array.display.matplotlib_qt import plt
from numina.array.display.matplotlib_qt import patches as patches
from numina.array.display.pause_debugplot import pause_debugplot
from emirdrp.instrument.csu_configuration import CsuConfiguration

from emirdrp.core import EMIR_NBARS
from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def display_slitlet_arrangement(fileobj,
                                bbox=None,
                                adjust=None,
                                geometry=None,
                                debugplot=0):
    """Display slitlet arrangment from CSUP keywords in FITS header.

    Parameters
    ----------
    fileobj : file object
        FITS file object.
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

    # read input FITS file
    hdulist = fits.open(fileobj.name)
    image_header = hdulist[0].header
    hdulist.close()

    # additional info from header
    grism = image_header['grism']
    spfilter = image_header['filter']
    rotang = image_header['rotang']

    # define slitlet arrangement
    csu_config = CsuConfiguration.define_from_fits(fileobj)

    # display arrangement
    if debugplot >= 10:
        print("slit     left    right   center   width")
        print("====  =======  =======  =======   =====")
        for i in range(EMIR_NBARS):
            ibar = i + 1
            print("{0:4d} {1:8.3f} {2:8.3f} {3:8.3f} {4:7.3f}".format(
                ibar, csu_config.csu_bar_left(ibar),
                csu_config.csu_bar_right(ibar),
                csu_config.csu_bar_slit_center(ibar),
                csu_config.csu_bar_slit_width(ibar)))
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
    if debugplot % 10 != 0:
        fig = plt.figure()
        if geometry is not None:
            x_geom, y_geom, dx_geom, dy_geom = geometry
            mngr = plt.get_current_fig_manager()
            mngr.window.setGeometry(x_geom, y_geom, dx_geom, dy_geom)
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
                ax.set_xlim([xmin, xmax])
            else:
                ax.set_xlim([0., 341.5])
            ax.set_ylim([0, 56])
        else:
            ax.set_xlim([bbox[0], bbox[1]])
            ax.set_ylim([bbox[2], bbox[3]])
        ax.set_xlabel('csu_bar_position (mm)')
        ax.set_ylabel('slit number')
        for i in range(EMIR_NBARS):
            ibar = i + 1
            ax.add_patch(patches.Rectangle(
                (csu_config.csu_bar_left(ibar), ibar-0.5),
                csu_config.csu_bar_slit_width(ibar), 1.0))
            ax.plot([0., csu_config.csu_bar_left(ibar)], [ibar, ibar], 'o-')
            ax.plot([csu_config.csu_bar_right(ibar), 341.5],
                    [ibar, ibar], 'o-')
        plt.title("File: " + fileobj.name + "\ngrism=" + grism +
                  ", filter=" + spfilter + ", rotang=" + str(round(rotang, 2)))
        pause_debugplot(debugplot, pltshow=True)

    # return results
    return csu_config._csu_bar_left, csu_config._csu_bar_right, \
           csu_config._csu_bar_slit_center, csu_config._csu_bar_slit_width


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("filename",
                        help="FITS files (wildcards accepted) or single TXT "
                             "file with list of FITS files",
                        type=argparse.FileType('r'),
                        nargs='+')

    # optional arguments
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

    list_fits_file_objects = []
    # if input file is a single txt file, assume it is a list of FITS files
    if len(args.filename) == 1:
        if args.filename[0].name[-4:] == ".txt":
            file_content = args.filename[0].read().splitlines()
            for line in file_content:
                if len(line) > 0:
                    if line[0] != '#':
                        tmpfile = line.split()[0]
                        if not os.path.isfile(tmpfile):
                            raise ValueError("File " + tmpfile + " not found!")
                        list_fits_file_objects.append(open(tmpfile, 'r'))
        else:
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
            display_slitlet_arrangement(fileobj, bbox=bbox,
                                        adjust=args.adjust,
                                        geometry=geometry,
                                        debugplot=args.debugplot)

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
