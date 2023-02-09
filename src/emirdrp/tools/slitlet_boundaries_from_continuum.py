#
# Copyright 2008-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
from datetime import datetime
import json
import numpy as np
import os.path
from scipy import ndimage
from scipy.signal import savgol_filter
import socket
import sys
from uuid import uuid4

from numina.array.ccd_line import SpectrumTrail
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow

from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.instrument.dtuconf import DtuConf

from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_VALID_GRISMS
from emirdrp.core import EMIR_VALID_FILTERS
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


class SlitletLimits(object):
    """Relevant parameters to fit boundaries of a given slitlet.

    Parameters
    ----------
    grism : string
        Grism name. It must be one in EMIR_VALID_GRISMS.
    spfilter : string
        Filter name. It must be one in EMIR_VALID_FILTERS.
    islitlet : int
        Number of slitlet.

    Attributes
    ----------
    grism : string
        Grism name. It must be one in EMIR_VALID_GRISMS.
    spfilter : string
        Filter name. It must be one in EMIR_VALID_FILTERS.
    islitlet : int
        Number of slitlet.
    bb_nc1_orig : int
        Minimum X coordinate (in pixel units) of the rectangle enclosing
        the 2D slitlet image.
    bb_nc2_orig : int
        Maximum X coordinate (in pixel units) of the rectangle enclosing
        the 2D slitlet image.
    bb_ns1_orig : int
        Minimum Y coordinate (in pixel units) of the rectangle enclosing
        the 2D slitlet image.
    bb_ns2_orig : int
        Maximum Y coordinate (in pixel units) of the rectangle enclosing
        the 2D slitlet image.
    xmin_lower_boundary_fit : int
        Minimum abscissa to fit lower boundary.
    xmax_lower_boundary_fit : int
        Maximum abscissa to fit lower boundary.
    xmin_upper_boundary_fit : int
        Minimum abscissa to fit upper boundary.
    xmax_upper_boundary_fit : int
        Maximum abscissa to fit upper boundary.

    """

    def __init__(self, grism, spfilter, islitlet):
        # parameters
        self.grism = grism
        self.spfilter = spfilter
        self.islitlet = islitlet

        # general protections
        if grism not in EMIR_VALID_GRISMS:
            raise ValueError("Grism " + str(grism) + " is not a valid option")
        if spfilter not in EMIR_VALID_FILTERS:
            raise ValueError("Filter " + str(spfilter) +
                             " is not a valid  option")
        if islitlet < 1 or islitlet > EMIR_NBARS:
            raise ValueError("islitlet=" + str(islitlet) +
                             " is outside valid range")

        # check that islitlet is within valid the range for the
        # grism + filter configuration
        if grism == "J" and spfilter == "J":
            if islitlet < 2 or islitlet > 54:
                raise ValueError("islitlet=" + str(islitlet) +
                                 " is outside valid range for grism " +
                                 str(grism) + " and filter " + str(spfilter))
        elif grism == "H" and spfilter == "H":
            if islitlet < 2 or islitlet > 54:
                raise ValueError("islitlet=" + str(islitlet) +
                                 " is outside valid range for grism " +
                                 str(grism) + " and filter " + str(spfilter))
        elif grism == "K" and spfilter == "Ksp":
            if islitlet < 2 or islitlet > 54:
                raise ValueError("islitlet=" + str(islitlet) +
                                 " is outside valid range for grism " +
                                 str(grism) + " and filter " + str(spfilter))
        elif grism == "LR" and spfilter == "YJ":
            if islitlet < 4 or islitlet > 55:
                raise ValueError("islitlet=" + str(islitlet) +
                                 " is outside valid range for grism " +
                                 str(grism) + " and filter " + str(spfilter))
        elif grism == "LR" and spfilter == "HK":
            if islitlet < 4 or islitlet > 55:
                raise ValueError("islitlet=" + str(islitlet) +
                                 " is outside valid range for grism " +
                                 str(grism) + " and filter " + str(spfilter))
        else:
            raise ValueError("Minimum and maximum islitlet still undefined "
                             "for grism " + str(grism) +
                             " and filter " + str(spfilter))

        # expected boundaries of the rectangle enclosing the 2D image
        coeff_bb_ns1 =  [-8.03677111e+01, 3.98169266e+01, -7.77949391e-02,
                         9.00823598e-04]
        delta_bb_ns2 = 84
        # offset measured overplotting ds9_bounddict.reg (created with the
        # script fit_boundaries) corresponding to grism J and filter J
        # over the same slitlet configuration with the corresponding grism
        # and filter setting
        if grism == "J" and spfilter == "J":
            offset_with_J_J = 0
        elif grism == "H" and spfilter == "H":
            offset_with_J_J = 3.0
        elif grism == "K" and spfilter == "Ksp":
            offset_with_J_J = 3.0
        elif grism == "LR" and spfilter == "YJ":
            offset_with_J_J = -95.0
        elif grism == "LR" and spfilter == "HK":
            offset_with_J_J = -95.0
        else:
            raise ValueError("Boundaries still undefined for grism " +
                             str(grism) + " and filter " + str(spfilter))
        coeff_bb_ns1[0] += offset_with_J_J
        poly_bb_ns1 = np.polynomial.Polynomial(coeff_bb_ns1)
        self.bb_nc1_orig = 1
        self.bb_nc2_orig = EMIR_NAXIS1

        self.bb_ns1_orig = int(poly_bb_ns1(islitlet) + 0.5)
        if self.bb_ns1_orig < 1:
            self.bb_ns1_orig = 1
        if self.bb_ns1_orig > EMIR_NAXIS2:
            raise ValueError("Unexpected bb_ns1_orig > EMIR_NAXIS2")
        self.bb_ns2_orig = int(poly_bb_ns1(islitlet) + 0.5) + delta_bb_ns2
        if self.bb_ns2_orig > EMIR_NAXIS2:
            self.bb_ns2_orig = EMIR_NAXIS2
        if self.bb_ns2_orig < 1:
            raise ValueError("Unexpected bb_ns2_orig < 1")

        # X ranges to fit boundaries
        self.xmin_lower_boundary_fit = 1
        self.xmax_lower_boundary_fit = EMIR_NAXIS1
        self.xmin_upper_boundary_fit = 1
        self.xmax_upper_boundary_fit = EMIR_NAXIS1

        # specific corrections for each grism + filter configuration
        if grism == "J" and spfilter == "J":
            self.deg_boundary = 5
            if islitlet == 54:
                self.xmin_upper_boundary_fit = 400
                self.xmax_upper_boundary_fit = 1750
        elif grism == "H" and spfilter == "H":
            self.deg_boundary = 5
            if islitlet == 54:
                self.xmin_upper_boundary_fit = 600
                self.xmax_upper_boundary_fit = 1600
        elif grism == "K" and spfilter == "Ksp":
            self.deg_boundary = 5
            if islitlet == 54:
                self.xmin_upper_boundary_fit = 400
                self.xmax_upper_boundary_fit = 1750
        elif grism == "LR" and spfilter == "YJ":
            self.deg_boundary = 5
            if islitlet == 4:
                self.xmin_lower_boundary_fit = 300
                # self.xmax_lower_boundary_fit = 1750
        elif grism == "LR" and spfilter == "HK":
            self.deg_boundary = 5
            if islitlet == 4:
                self.xmin_lower_boundary_fit = 300
                # self.xmax_lower_boundary_fit = 1750
        else:
            raise ValueError("Ranges to fit boundaries still undefined "
                             "for grism " + str(grism) +
                             " and filter " + str(spfilter))


def extract_slitlet2d(image_2k2k, sltlim):
    """Extract slitlet 2d image from image with original EMIR dimensions.

    Parameters
    ----------
    image_2k2k : 2d numpy array, float
        Original image (dimensions NAXIS1 * NAXIS2)
    sltlim : instance of SlitLimits class
        Object containing relevant information concerning the slitlet
        region to be extracted.

    Returns
    -------
    slitlet2d : 2d numpy array, float
        Image corresponding to the slitlet region defined by its
        bounding box.

    """

    # extract slitlet region
    slitlet2d = image_2k2k[(sltlim.bb_ns1_orig - 1):sltlim.bb_ns2_orig,
                           (sltlim.bb_nc1_orig - 1):sltlim.bb_nc2_orig]

    # transform to float
    slitlet2d = slitlet2d.astype(float)

    # return slitlet image
    return slitlet2d


def compute_slitlet_boundaries(
        filename, grism, spfilter, list_slitlets,
        size_x_medfilt, size_y_savgol,
        times_sigma_threshold,
        bounddict, debugplot=0):
    """Compute slitlet boundaries using continuum lamp images.

    Parameters
    ----------
    filename : string
        Input continuum lamp image.
    grism : string
        Grism name. It must be one in EMIR_VALID_GRISMS.
    spfilter : string
        Filter name. It must be one in EMIR_VALID_FILTERS.
    list_slitlets : list of integers
        Number of slitlets to be updated.
    size_x_medfilt : int
        Window in the X (spectral) direction, in pixels, to apply
        the 1d median filter in order to remove bad pixels.
    size_y_savgol : int
        Window in the Y (spatial) direction to be used when using the
        1d Savitzky-Golay filter.
    times_sigma_threshold : float
        Times sigma to detect peaks in derivatives.
    bounddict : dictionary of dictionaries
        Structure to store the boundaries.
    debugplot : int
        Determines whether intermediate computations and/or plots are
        displayed.

    """

    # read 2D image
    hdulist = fits.open(filename)
    image_header = hdulist[0].header
    image2d = hdulist[0].data
    naxis2, naxis1 = image2d.shape
    hdulist.close()
    if debugplot >= 10:
        print('>>> NAXIS1:', naxis1)
        print('>>> NAXIS2:', naxis2)


    # ToDo: replace this by application of cosmetic defect mask!
    for j in range(1024):
        image2d[1024, j] = (image2d[1023, j] + image2d[1025, j]) / 2
        image2d[1023, j + 1024] = (image2d[1022, j + 1024] +
                                   image2d[1024, j + 1024]) / 2
        
    # remove path from filename
    sfilename = os.path.basename(filename)

    # check that the FITS file has been obtained with EMIR
    instrument = image_header['instrume']
    if instrument != 'EMIR':
        raise ValueError("INSTRUME keyword is not 'EMIR'!")

    # read CSU configuration from FITS header
    csu_config = CsuConfiguration.define_from_fits(filename)
    # read DTU configuration from FITS header
    with fits.open(filename) as img:
        dtu2_config = DtuConf.from_img(img)

    # read grism
    grism_in_header = image_header['grism']
    if grism != grism_in_header:
        raise ValueError("GRISM keyword=" + grism_in_header +
                         " is not the expected value=" + grism)

    # read filter
    spfilter_in_header = image_header['filter']
    if spfilter != spfilter_in_header:
        raise ValueError("FILTER keyword=" + spfilter_in_header +
                         " is not the expected value=" + spfilter)

    # read rotator position angle
    rotang = image_header['rotang']

    # read date-obs
    date_obs = image_header['date-obs']

    for islitlet in list_slitlets:
        if debugplot < 10:
            sys.stdout.write('.')
            sys.stdout.flush()

        sltlim = SlitletLimits(grism, spfilter, islitlet)
        # extract slitlet2d
        slitlet2d = extract_slitlet2d(image2d, sltlim)
        if debugplot % 10 != 0:
            ximshow(slitlet2d,
                    title=sfilename + " [original]"
                          "\nslitlet=" + str(islitlet) +
                          ", grism=" + grism +
                          ", filter=" + spfilter +
                          ", rotang=" + str(round(rotang, 2)),
                    first_pixel=(sltlim.bb_nc1_orig, sltlim.bb_ns1_orig),
                    debugplot=debugplot)

        # apply 1d median filtering (along the spectral direction)
        # to remove bad pixels
        size_x = size_x_medfilt
        size_y = 1
        slitlet2d_smooth = ndimage.filters.median_filter(
            slitlet2d, size=(size_y, size_x))

        if debugplot % 10 != 0:
            ximshow(slitlet2d_smooth,
                    title=sfilename + " [smoothed]"
                          "\nslitlet=" + str(islitlet) +
                          ", grism=" + grism +
                          ", filter=" + spfilter +
                          ", rotang=" + str(round(rotang, 2)),
                    first_pixel=(sltlim.bb_nc1_orig, sltlim.bb_ns1_orig),
                    debugplot=debugplot)

        # apply 1d Savitzky-Golay filter (along the spatial direction)
        # to compute first derivative
        slitlet2d_savgol = savgol_filter(
            slitlet2d_smooth, window_length=size_y_savgol, polyorder=2,
            deriv=1, axis=0)

        # compute basic statistics
        q25, q50, q75 = np.percentile(slitlet2d_savgol, q=[25.0, 50.0, 75.0])
        sigmag = 0.7413 * (q75 - q25)  # robust standard deviation
        if debugplot >= 10:
            print("q50, sigmag:", q50, sigmag)

        if debugplot % 10 != 0:
            ximshow(slitlet2d_savgol,
                    title=sfilename + " [S.-G.filt.]"
                          "\nslitlet=" + str(islitlet) +
                          ", grism=" + grism +
                          ", filter=" + spfilter +
                          ", rotang=" + str(round(rotang, 2)),
                    first_pixel=(sltlim.bb_nc1_orig, sltlim.bb_ns1_orig),
                    z1z2=(q50-times_sigma_threshold*sigmag,
                          q50+times_sigma_threshold*sigmag),
                    debugplot=debugplot)

        # identify objects in slitlet2d_savgol: pixels with positive
        # derivatives are identify independently from pixels with
        # negative derivaties; then the two set of potential features
        # are merged; this approach avoids some problems when, in
        # nearby regions, there are pixels with positive and negative
        # derivatives (in those circumstances a single search as
        # np.logical_or(
        #     slitlet2d_savgol < q50 - times_sigma_threshold * sigmag,
        #     slitlet2d_savgol > q50 + times_sigma_threshold * sigmag)
        # led to erroneous detections!)
        #
        # search for positive derivatives
        labels2d_objects_pos, no_objects_pos = ndimage.label(
            slitlet2d_savgol > q50 + times_sigma_threshold * sigmag)
        # search for negative derivatives
        labels2d_objects_neg, no_objects_neg = ndimage.label(
            slitlet2d_savgol < q50 - times_sigma_threshold * sigmag)
        # merge both sets
        non_zero_neg = np.where(labels2d_objects_neg > 0)
        labels2d_objects = np.copy(labels2d_objects_pos)
        labels2d_objects[non_zero_neg] += \
            labels2d_objects_neg[non_zero_neg] + no_objects_pos
        no_objects = no_objects_pos + no_objects_neg

        if debugplot >= 10:
            print("Number of objects with positive derivative:",
                  no_objects_pos)
            print("Number of objects with negative derivative:",
                  no_objects_neg)
            print("Total number of objects initially found...:", no_objects)

        if debugplot % 10 != 0:
            ximshow(labels2d_objects,
                    z1z2=(0, no_objects),
                    title=sfilename + " [objects]"
                                     "\nslitlet=" + str(islitlet) +
                          ", grism=" + grism +
                          ", filter=" + spfilter +
                          ", rotang=" + str(round(rotang, 2)),
                    first_pixel=(sltlim.bb_nc1_orig, sltlim.bb_ns1_orig),
                    cbar_label="Object number",
                    debugplot=debugplot)

        # select boundaries as the largest objects found with
        # positive and negative derivatives
        n_der_pos = 0  # number of pixels covered by the object with deriv > 0
        i_der_pos = 0  # id of the object with deriv > 0
        n_der_neg = 0  # number of pixels covered by the object with deriv < 0
        i_der_neg = 0  # id of the object with deriv < 0
        for i in range(1, no_objects+1):
            xy_tmp = np.where(labels2d_objects == i)
            n_pix = len(xy_tmp[0])
            if i <= no_objects_pos:
                if n_pix > n_der_pos:
                    i_der_pos = i
                    n_der_pos = n_pix
            else:
                if n_pix > n_der_neg:
                    i_der_neg = i
                    n_der_neg = n_pix

        # determine which boundary is lower and which is upper
        y_center_mass_der_pos = ndimage.center_of_mass(
            slitlet2d_savgol, labels2d_objects, [i_der_pos])[0][0]
        y_center_mass_der_neg = ndimage.center_of_mass(
            slitlet2d_savgol, labels2d_objects, [i_der_neg])[0][0]
        if y_center_mass_der_pos < y_center_mass_der_neg:
            i_lower = i_der_pos
            i_upper = i_der_neg
            if debugplot >= 10:
                print("-> lower boundary has positive derivatives")
        else:
            i_lower = i_der_neg
            i_upper = i_der_pos
            if debugplot >= 10:
                print("-> lower boundary has negative derivatives")
        list_slices_ok = [i_lower, i_upper]

        # adjust individual boundaries passing the selection:
        # - select points in the image belonging to a given boundary
        # - compute weighted mean of the pixels of the boundary, column
        #   by column (this reduces dramatically the number of points
        #   to be fitted to determine the boundary)
        list_boundaries = []
        for k in range(2):  # k=0 lower boundary, k=1 upper boundary
            # select points to be fitted for a particular boundary
            # (note: be careful with array indices and pixel
            # coordinates)
            xy_tmp = np.where(labels2d_objects == list_slices_ok[k])
            xmin = xy_tmp[1].min()  # array indices (integers)
            xmax = xy_tmp[1].max()  # array indices (integers)
            xfit = []
            yfit = []
            # fix range for fit
            if k == 0:
                xmineff = max(sltlim.xmin_lower_boundary_fit, xmin)
                xmaxeff = min(sltlim.xmax_lower_boundary_fit, xmax)
            else:
                xmineff = max(sltlim.xmin_upper_boundary_fit, xmin)
                xmaxeff = min(sltlim.xmax_upper_boundary_fit, xmax)
            # loop in columns of the image belonging to the boundary
            for xdum in range(xmineff, xmaxeff + 1):  # array indices (integer)
                iok = np.where(xy_tmp[1] == xdum)
                y_tmp = xy_tmp[0][iok] + sltlim.bb_ns1_orig  # image pixel
                weight = slitlet2d_savgol[xy_tmp[0][iok], xy_tmp[1][iok]]
                y_wmean = sum(y_tmp * weight) / sum(weight)
                xfit.append(xdum + sltlim.bb_nc1_orig)
                yfit.append(y_wmean)
            xfit = np.array(xfit)
            yfit = np.array(yfit)
            # declare new SpectrumTrail instance
            boundary = SpectrumTrail()
            # define new boundary
            boundary.fit(x=xfit, y=yfit, deg=sltlim.deg_boundary,
                         times_sigma_reject=10,
                         title="slit:" + str(sltlim.islitlet) +
                               ", deg=" + str(sltlim.deg_boundary),
                         debugplot=0)
            list_boundaries.append(boundary)

        if debugplot % 10 != 0:
            for tmp_img, tmp_label in zip(
                [slitlet2d_savgol, slitlet2d],
                [' [S.-G.filt.]', ' [original]']
            ):
                ax = ximshow(tmp_img,
                             title=sfilename + tmp_label +
                                   "\nslitlet=" + str(islitlet) +
                                   ", grism=" + grism +
                                   ", filter=" + spfilter +
                                   ", rotang=" + str(round(rotang, 2)),
                             first_pixel=(sltlim.bb_nc1_orig,
                                          sltlim.bb_ns1_orig),
                             show=False,
                             debugplot=debugplot)
                for k in range(2):
                    xpol, ypol = list_boundaries[k].linspace_pix(
                        start=1, stop=EMIR_NAXIS1)
                    ax.plot(xpol, ypol, 'b--', linewidth=1)
                for k in range(2):
                    xpol, ypol = list_boundaries[k].linspace_pix()
                    ax.plot(xpol, ypol, 'g--', linewidth=4)
                # show plot
                pause_debugplot(debugplot, pltshow=True)

        # update bounddict
        tmp_dict = {
            'boundary_coef_lower':
                list_boundaries[0].poly_funct.coef.tolist(),
            'boundary_xmin_lower': list_boundaries[0].xlower_line,
            'boundary_xmax_lower': list_boundaries[0].xupper_line,
            'boundary_coef_upper':
                list_boundaries[1].poly_funct.coef.tolist(),
            'boundary_xmin_upper': list_boundaries[1].xlower_line,
            'boundary_xmax_upper': list_boundaries[1].xupper_line,
            'csu_bar_left': csu_config.csu_bar_left(islitlet),
            'csu_bar_right': csu_config.csu_bar_right(islitlet),
            'csu_bar_slit_center': csu_config.csu_bar_slit_center(islitlet),
            'csu_bar_slit_width': csu_config.csu_bar_slit_width(islitlet),
            'rotang': rotang,
            'xdtu': dtu2_config.xaxis.coor,
            'ydtu': dtu2_config.yaxis.coor,
            'zdtu': dtu2_config.zaxis.coor,
            'xdtu_0': dtu2_config.xaxis.coorx_0,
            'ydtu_0': dtu2_config.yaxis.coory_0,
            'zdtu_0': dtu2_config.zaxis.coor_0,
            'zzz_info1': os.getlogin() + '@' + socket.gethostname(),
            'zzz_info2': datetime.now().isoformat()
        }
        slitlet_label = "slitlet" + str(islitlet).zfill(2)
        if slitlet_label not in bounddict['contents']:
            bounddict['contents'][slitlet_label] = {}
        bounddict['contents'][slitlet_label][date_obs] = tmp_dict

    if debugplot < 10:
        print("")


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("filename",
                        help="FITS file or txt file with list of FITS files",
                        type=argparse.FileType('rb'))
    parser.add_argument("--grism", required=True,
                        help="Grism name",
                        choices=EMIR_VALID_GRISMS)
    parser.add_argument("--filter", required=True,
                        help="Filter name",
                        choices=EMIR_VALID_FILTERS)
    parser.add_argument("--tuple_slit_numbers", required=True,
                        help="Tuple n1[,n2[,step]] to define slitlet numbers")

    # optional arguments
    parser.add_argument("--first_time",
                        help="Generate new bounddict json file",
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

    # read slitlet numbers to be computed
    tmp_str = args.tuple_slit_numbers.split(",")
    if len(tmp_str) == 3:
        n1 = int(tmp_str[0])
        n2 = int(tmp_str[1])
        step = int(tmp_str[2])
    elif len(tmp_str) == 2:
        n1 = int(tmp_str[0])
        n2 = int(tmp_str[1])
        step = 1
    elif len(tmp_str) == 1:
        n1 = int(tmp_str[0])
        n2 = n1
        step = 1
    else:
        raise ValueError("Invalid tuple for slitlet numbers")
    if n1 < 1:
        raise ValueError("Invalid slitlet number < 1")
    if n2 > EMIR_NBARS:
        raise ValueError("Invalid slitlet number > EMIR_NBARS")
    if step <= 0:
        raise ValueError("Invalid step <= 0")
    list_slitlets = range(n1, n2+1, step)

    # define bounddict file name
    bounddict_file = "bounddict_grism_" + args.grism + "_filter_" + \
                     args.filter + ".json"

    # define bounddict prior to the new computation
    if args.first_time:
        bounddict = {}
        bounddict['instrument'] = 'EMIR'
        bounddict['meta_info'] = {}
        bounddict['meta_info']['creation_date'] = datetime.now().isoformat()
        bounddict['meta_info']['description'] = \
            'slitlet boundaries from fits to continuum-lamp exposures'
        bounddict['meta_info']['recipe_name'] = 'undefined'
        bounddict['tags'] = {}
        bounddict['tags']['grism'] = args.grism
        bounddict['tags']['filter'] = args.filter
        bounddict['uuid'] = str(uuid4())
        bounddict['contents'] = {}
        print('>>> Generating new bounddict from scratch')
    else:
        if not os.path.isfile(bounddict_file):
            raise ValueError("File " + bounddict_file + " not found!")
        else:
            bounddict = json.loads(open(bounddict_file, mode='rt').read())
            print('>>> Initializing bounddict from previous file:')
            print(bounddict_file)

    # if input file is a txt file, assume it is a list of FITS files
    filename = args.filename.name
    if filename[-4:] == ".txt":
        # since the input filename has been opened with argparse in binary
        # mode, it is necessary to close it and open it in text mode
        args.filename.close()
        with open(filename, mode='rt') as f:
            file_content = f.read().splitlines()
        list_fits_files = []
        for line in file_content:
            if len(line) > 0:
                if line[0] != '#':
                    tmpfile = line.split()[0]
                    if not os.path.isfile(tmpfile):
                        raise ValueError("File " + tmpfile + " not found!")
                    list_fits_files.append(tmpfile)
    else:
        list_fits_files = [filename]

    # update bounddict
    for ifile, myfile in enumerate(list_fits_files):
        print('>>> Reading file ' + str(ifile+1) + "/" +
              str(len(list_fits_files)) + ":\n" + myfile)
        compute_slitlet_boundaries(
            filename=myfile,
            list_slitlets=list_slitlets,
            grism=args.grism, spfilter=args.filter,
            size_x_medfilt=31,
            size_y_savgol=11,
            times_sigma_threshold=3.0,
            bounddict=bounddict, debugplot=args.debugplot)

    # save new version of bounddict
    # print(json.dumps(bounddict, indent=4, sort_keys=True))
    with open(bounddict_file, 'wt') as fstream:
        json.dump(bounddict, fstream, indent=2, sort_keys=True)
        print('>>> Saving file ' + bounddict_file)


if __name__ == "__main__":

    main()
