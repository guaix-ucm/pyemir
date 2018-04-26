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


"""Wavelength calibration functionality"""

from astropy.io import fits
from datetime import datetime
import logging
from matplotlib.patches import Rectangle
import numpy as np
from scipy import ndimage
from scipy.interpolate import interp1d
from scipy.signal import medfilt
from skimage import restoration
import sys
import time
from uuid import uuid4

from numina.array.ccd_line import ArcLine
from numina.array.ccd_line import intersection_spectrail_arcline
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.polfit_residuals import \
    polfit_residuals_with_sigma_rejection
from numina.array.display.ximplotxy import ximplotxy
from numina.array.display.ximshow import ximshow
from numina.array.distortion import compute_distortion
from numina.array.distortion import fmap
from numina.array.distortion import ncoef_fmap
from numina.array.distortion import order_fmap
from numina.array.distortion import rectify2d
from numina.array.wavecalib.__main__ import read_wv_master_from_array
from numina.array.wavecalib.__main__ import wvcal_spectrum
from numina.array.wavecalib.arccalibration import refine_arccalibration
from numina.array.wavecalib.peaks_spectrum import find_peaks_spectrum
from numina.array.wavecalib.peaks_spectrum import refine_peaks_spectrum
from numina.array.wavecalib.resample import resample_image2d_flux
import numina.types.qc

from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.instrument.dtu_configuration import DtuConfiguration
from emirdrp.products import RectWaveCoeff
from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.fit_boundaries import expected_distorted_boundaries
from emirdrp.tools.fit_boundaries import expected_distorted_frontiers
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers
from emirdrp.tools.select_unrectified_slitlets import \
    select_unrectified_slitlet

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_VALID_FILTERS
from emirdrp.core import EMIR_VALID_GRISMS


class Slitlet2dArc(object):
    """Slitlet2dArc class definition.

    Slitlet2dArc: 2D slitlet class for Long-Slit Arc image

    It is important to distinguish between boundaries (the slitlet
    region when useful information is available) and frontiers (which
    define the real separation between consecutive slitlets when no gap
    between them is considered).

    Parameters
    ----------
    islitlet : int
        Slitlet number.
    params : :class:`~lmfit.parameter.Parameters`
        Parameters to be employed in the prediction of the distorted
        boundaries.
    parmodel : str
        Model to be assumed. Allowed values are 'longslit' and
        'multislit'.
    csu_conf : CsuConfiguration object
        Instance of CsuConfiguration.
    ymargin_bb : int
        Extra number of pixels above and below the enclosing rectangle
        defined by the slitlet frontiers.
    debugplot : int
        Debugging level for messages and plots. For details see
        'numina.array.display.pause_debugplot.py'.

   Attributes
    ----------
    islitlet : int
        Slitlet number.
    csu_bar_left : float
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
    ymargin_bb : int
        Extra number of pixels above and below the enclosing rectangle
        defined by the slitlet frontiers.
    bb_nc1_orig : int
        Minimum X coordinate of the enclosing bounding box (in pixel
        units) in the original image.
    bb_nc2_orig : int
        Maximum X coordinate of the enclosing bounding box (in pixel
        units) in the original image.
    bb_ns1_orig : int
        Minimum Y coordinate of the enclosing bounding box (in pixel
        units) in the original image.
    bb_ns2_orig : int
        Maximum Y coordinate of the enclosing bounding box (in pixel
        units) in the original image.
    x0_reference : float
        X coordinate where the rectified y0_reference_middle is computed
        as the Y coordinate of the middle spectrum trail. The same value
        is used for all the available spectrum trails.
    list_spectrails: list of SpectrumTrail instances
        List of spectrum trails defined.
    list_frontiers: list of SpectrumTrail instances
        List of spectrum trails defining the slitlet frontiers.
    y0_reference_lower: float
        Y coordinate corresponding to the lower spectrum trail computed
        at x0_reference. This value is employed as the Y coordinate of
        the lower spectrum trail of the rectified slitlet.
    y0_reference_middle: float
        Y coordinate corresponding to the middle spectrum trail computed
        at x0_reference. This value is employed as the Y coordinate of
        the middle spectrum trail of the rectified slitlet.
    y0_reference_upper: float
        Y coordinate corresponding to the upper spectrum trail computed
        at x0_reference. This value is employed as the Y coordinate of
        the upper spectrum trail of the rectified slitlet.
    y0_frontier_lower: float
        Y coordinate corresponding to the lower frontier computed at
        x0_reference.
    y0_frontier_upper: float
        Y coordinate corresponding to the upper frontier computed at
        x0_reference.
    x_inter_orig : 1d numpy array, float
        X coordinates of the intersection points of arc lines with
        spectrum trails in the original image.
    y_inter_orig : 1d numpy array, float
        Y coordinates of the intersection points of arc lines with
        spectrum trails in the original image.
    x_inter_rect : 1d numpy array, float
        X coordinates of the intersection points of arc lines with
        spectrum trails in the rectified image.
    y_inter_rect : 1d numpy array, float
        Y coordinates of the intersection points of arc lines with
        spectrum trails in the rectified image.
    ttd_order : int or None
        Polynomial order corresponding to the rectification
        transformation. It is None until the coefficients ttd_aij and
        ttd_bij have been computed.
    ttd_aij : numpy array
        Polynomial coefficents corresponding to the direct
        rectification transformation coefficients a_ij.
    ttd_bij : numpy array
        Polynomial coefficents corresponding to the direct
        rectification transformation coefficients b_ij.
    tti_aij : numpy array
        Polynomial coefficents corresponding to the inverse
        rectification transformation coefficients a_ij.
    tti_bij : numpy array
        Polynomial coefficents corresponding to the inverse
        rectification transformation coefficients b_ij.
    ttd_order_longslit_model : int or None
        Polynomial order corresponding to the modeled direct
        rectification transformation. It is None until the coefficients
        ttd_aij_longslit_model and ttd_bij_longslit_model have been
        computed.
    ttd_aij_longslit_model : numpy array
        Polynomial coefficents corresponding to the direct
        rectification transformation coefficients a_ij interpolated
        with a smooth polynomial variation as a function of
        y0_reference_middle.
    ttd_bij_longslit_model : numpy array
        Polynomial coefficents corresponding to the direct
        rectification transformation coefficients b_ij interpolated
        with a smooth polynomial variation as a function of
        y0_reference_middle.
    tti_aij_longslit_model : numpy array
        Polynomial coefficents corresponding to the inverse
        rectification transformation coefficients a_ij interpolated
        with a smooth polynomial variation as a function of
        y0_reference_middle.
    tti_bij_longslit_model : numpy array
        Polynomial coefficents corresponding to the inverse
        rectification transformation coefficients b_ij interpolated
        with a smooth polynomial variation as a function of
        y0_reference_middle.
    wpoly : numpy.polynomial.Polynomial instance
        Wavelength calibration polynomial, providing the
        wavelength as a function of pixel number (running from 1 to
        NAXIS1).
    wpoly_longslit_model : numpy.polynomial.Polynomial instance
        Refined and modeled wavelength calibration polynomial,
        providing the wavelength as a function of pixel number (running
        from 1 to NAXIS1), or None (when the fit cannot be obtained).
    crval1_linear : float
        CRVAL1 corresponding to a linear wavelength scale computed from
        wpoly.
    cdelt1_linear : float
        CDELT1 corresponding to a linear wavelength scale computed from
        wpoly.
    debugplot : int
        Debugging level for messages and plots. For details see
        'numina.array.display.pause_debugplot.py'.

    """

    def __init__(self, islitlet, params, parmodel, csu_conf, ymargin_bb,
                 debugplot):

        # slitlet number
        self.islitlet = islitlet

        # csu configuration
        self.csu_bar_left = csu_conf.csu_bar_left(islitlet)
        self.csu_bar_right = csu_conf.csu_bar_right(islitlet)
        self.csu_bar_slit_center = csu_conf.csu_bar_slit_center(islitlet)
        self.csu_bar_slit_width = csu_conf.csu_bar_slit_width(islitlet)

        # horizontal bounding box
        self.bb_nc1_orig = 1
        self.bb_nc2_orig = EMIR_NAXIS1

        # reference abscissa
        self.x0_reference = float(EMIR_NAXIS1) / 2.0 + 0.5  # single float

        # compute spectrum trails and store in which order they are computed
        self.i_lower_spectrail = 0
        self.i_middle_spectrail = 1
        self.i_upper_spectrail = 2
        self.list_spectrails = expected_distorted_boundaries(
                islitlet, self.csu_bar_slit_center,
                [0, 0.5, 1], params, parmodel,
                numpts=101, deg=5, debugplot=0
            )
        # update y_rectified computed at x0_reference
        for spectrail in self.list_spectrails:
            spectrail.y_rectified = spectrail.poly_funct(self.x0_reference)

        # define reference ordinates using lower, middle and upper spectrails
        # evaluated at x0_reference
        self.y0_reference_lower = \
            self.list_spectrails[self.i_lower_spectrail].y_rectified
        self.y0_reference_middle = \
            self.list_spectrails[self.i_middle_spectrail].y_rectified
        self.y0_reference_upper = \
            self.list_spectrails[self.i_upper_spectrail].y_rectified

        # compute frontiers (lower and upper)
        self.list_frontiers = expected_distorted_frontiers(
            islitlet, self.csu_bar_slit_center,
            params, parmodel,
            numpts=101, deg=5, debugplot=0
        )
        # update y_rectified computed at x0_reference
        for spectrail in self.list_frontiers:
            spectrail.y_rectified = spectrail.poly_funct(self.x0_reference)

        # define frontier ordinates at x0_reference
        self.y0_frontier_lower = self.list_frontiers[0].y_rectified
        self.y0_frontier_upper = self.list_frontiers[1].y_rectified

        # determine vertical bounding box
        self.ymargin_bb = ymargin_bb
        xdum = np.linspace(1, EMIR_NAXIS1, num=EMIR_NAXIS1)
        ylower = self.list_frontiers[0].poly_funct(xdum)
        yupper = self.list_frontiers[1].poly_funct(xdum)
        self.bb_ns1_orig = int(ylower.min() + 0.5) - self.ymargin_bb
        if self.bb_ns1_orig < 1:
            self.bb_ns1_orig = 1
        self.bb_ns2_orig = int(yupper.max() + 0.5) + self.ymargin_bb
        if self.bb_ns2_orig > EMIR_NAXIS2:
            self.bb_ns2_orig = EMIR_NAXIS2

        # place holder for still undefined members
        self.list_arc_lines = None
        self.x_inter_orig = None
        self.y_inter_orig = None
        self.x_inter_rect = None
        self.y_inter_rect = None
        self.ttd_order = None
        self.ttd_aij = None
        self.ttd_bij = None
        self.tti_aij = None
        self.tti_bij = None
        self.wpoly = None
        self.crval1_linear = None
        self.cdelt1_linear = None

        # undefined members to be updated with images in longslit mode
        self.ttd_order_longslit_model = None
        self.ttd_aij_longslit_model = None
        self.ttd_bij_longslit_model = None
        self.tti_aij_longslit_model = None
        self.tti_bij_longslit_model = None
        self.wpoly_longslit_model = None

        # debugplot
        self.debugplot = debugplot

    def __repr__(self):
        """Printable representation of a Slitlet2dArc instance."""

        # string with all the information
        output = "<Slitlet2dArc instance>\n" + \
            "- islitlet....................: " + \
                 str(self.islitlet) + "\n" + \
            "- csu_bar_left................: " + \
                 str(self.csu_bar_left) + "\n" + \
            "- csu_bar_right...............: " + \
                 str(self.csu_bar_right) + "\n" + \
            "- csu_bar_slit_center.........: " + \
                 str(self.csu_bar_slit_center) + "\n" + \
            "- csu_bar_slit_width..........: " + \
                 str(self.csu_bar_slit_width) + "\n" + \
            "- x0_reference................: " + \
                 str(self.x0_reference) + "\n" + \
            "- y0_reference_lower..........: " + \
                 str(self.y0_reference_lower) + "\n" + \
            "- y0_reference_middle.........: " + \
                 str(self.y0_reference_middle) + "\n" + \
            "- y0_reference_upper..........: " + \
                 str(self.y0_reference_upper) + "\n" + \
            "- y0_frontier_lower..........: " + \
                 str(self.y0_frontier_lower) + "\n" + \
            "- y0_frontier_upper..........: " + \
                 str(self.y0_frontier_upper) + "\n" + \
            "- bb_nc1_orig.................: " + \
                 str(self.bb_nc1_orig) + "\n" + \
            "- bb_nc2_orig.................: " + \
                 str(self.bb_nc2_orig) + "\n" + \
            "- bb_ns1_orig.................: " + \
                 str(self.bb_ns1_orig) + "\n" + \
            "- bb_ns2_orig.................: " + \
                 str(self.bb_ns2_orig) + "\n" + \
            "- lower spectrail_poly_funct..:\n\t" + \
                 str(self.list_spectrails[self.i_lower_spectrail].poly_funct)\
                 + "\n" + \
            "- middle spectrail_poly_funct.:\n\t" + \
                 str(self.list_spectrails[self.i_middle_spectrail].poly_funct)\
                 + "\n" + \
            "- upper spectrail_poly_funct..:\n\t" + \
                 str(self.list_spectrails[self.i_upper_spectrail].poly_funct)\
                 + "\n" + \
            "- lower frontier_poly_funct...:\n\t" + \
                 str(self.list_frontiers[0].poly_funct) + "\n" + \
            "- upper frontier_poly_funct...:\n\t" + \
                 str(self.list_frontiers[1].poly_funct) + "\n"

        if self.list_arc_lines is None:
            number_arc_lines = None
        else:
            number_arc_lines = len(self.list_arc_lines)

        output += "- num. of associated arc lines: " + \
                  str(number_arc_lines) + "\n"

        for dumval, dumlab in zip(
                [self.x_inter_orig, self.y_inter_orig,
                 self.x_inter_rect, self.y_inter_rect],
                ["x_inter_orig", "y_inter_orig",
                 "x_inter_rect", "y_inter_rect"]
        ):
            if dumval is None:
                output += "- " + dumlab + "................: None\n"
            else:
                output += "- " + dumlab + "................: " + \
                          str(len(dumval)) + " values defined\n\t[" + \
                          str(dumval[0]) + ", ..., " + \
                          str(dumval[-1]) + "]\n"
        output += \
            "- ttd_order...............: " + str(self.ttd_order) + "\n" + \
            "- ttd_aij.................:\n\t" + str(self.ttd_aij) + "\n" + \
            "- ttd_bij.................:\n\t" + str(self.ttd_bij) + "\n" + \
            "- tti_aij.................:\n\t" + str(self.tti_aij) + "\n" + \
            "- tti_bij.................:\n\t" + str(self.tti_bij) + "\n" + \
            "- ttd_order_longslit_model: " + \
            str(self.ttd_order_longslit_model) + "\n" + \
            "- ttd_aij_longslit_model..:\n\t" + \
            str(self.ttd_aij_longslit_model) + "\n" + \
            "- ttd_bij_longslit_model..:\n\t" + \
            str(self.ttd_bij_longslit_model) + "\n" + \
            "- tti_aij_longslit_model..:\n\t" + \
            str(self.tti_aij_longslit_model) + "\n" + \
            "- tti_bij_longslit_model..:\n\t" + \
            str(self.tti_bij_longslit_model) + "\n" + \
            "- wpoly..................:\n\t" + \
            str(self.wpoly) + "\n" + \
            "- wpoly_longslit_model....:\n\t" + \
            str(self.wpoly_longslit_model) + "\n" + \
            "- crval1_linear...........: " + str(self.crval1_linear) + "\n" + \
            "- cdelt1_linear...........: " + str(self.cdelt1_linear) + "\n" + \
            "- debugplot...............: " + \
            str(self.debugplot)

        return output

    def extract_slitlet2d(self, image_2k2k):
        """Extract slitlet 2d image from image with original EMIR dimensions.

        Parameters
        ----------
        image_2k2k : numpy array
            Original image (dimensions EMIR_NAXIS1 * EMIR_NAXIS2)

        Returns
        -------
        slitlet2d : numpy array
            Image corresponding to the slitlet region defined by its
            bounding box.

        """

        # protections
        naxis2, naxis1 = image_2k2k.shape
        if naxis1 != EMIR_NAXIS1:
            raise ValueError('Unexpected naxis1')
        if naxis2 != EMIR_NAXIS2:
            raise ValueError('Unexpected naxis2')

        # extract slitlet region
        slitlet2d = image_2k2k[(self.bb_ns1_orig - 1):self.bb_ns2_orig,
                               (self.bb_nc1_orig - 1):self.bb_nc2_orig]

        # transform to float
        slitlet2d = slitlet2d.astype(np.float)

        # display slitlet2d with boundaries and middle spectrum trail
        if abs(self.debugplot) in [21, 22]:
            ax = ximshow(slitlet2d, title="Slitlet#" + str(self.islitlet),
                         first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                         show=False)
            xdum = np.linspace(1, EMIR_NAXIS1, num=EMIR_NAXIS1)
            ylower = \
                self.list_spectrails[self.i_lower_spectrail].poly_funct(xdum)
            ax.plot(xdum, ylower, 'b-')
            ymiddle = \
                self.list_spectrails[self.i_middle_spectrail].poly_funct(xdum)
            ax.plot(xdum, ymiddle, 'b--')
            yupper = \
                self.list_spectrails[self.i_upper_spectrail].poly_funct(xdum)
            ax.plot(xdum, yupper, 'b-')
            ylower_frontier = self.list_frontiers[0].poly_funct(xdum)
            ax.plot(xdum, ylower_frontier, 'b:')
            yupper_frontier = self.list_frontiers[1].poly_funct(xdum)
            ax.plot(xdum, yupper_frontier, 'b:')
            pause_debugplot(debugplot=self.debugplot, pltshow=True)

        # return slitlet image
        return slitlet2d

    def locate_unknown_arc_lines(self, slitlet2d,
                                 times_sigma_threshold=15,
                                 minimum_threshold=None,
                                 delta_x_max=30,
                                 delta_y_min=30,
                                 min_dist_from_middle=15):
        """Determine the location of known arc lines in slitlet.

        Parameters
        ----------
        slitlet2d : numpy array
            Image containing the 2d slitlet image.
        times_sigma_threshold : float
            Times (robust) sigma above the median of the image to look
            for arc lines.
        minimum_threshold : float or None
            Minimum threshold to look for arc lines.
        delta_x_max : float
            Maximum size of potential arc line in the X direction.
        delta_y_min : float
            Minimum size of potential arc line in the Y direction.
        min_dist_from_middle : float
            Minimum Y distance from the middle spectrum trail to the
            extreme of the potential arc line. This constraint avoid
            detecting arc line reflections as bone fide arc lines.

        """

        # smooth denoising of slitlet2d
        slitlet2d_rs, coef_rs = rescale_array_to_z1z2(slitlet2d, z1z2=(-1, 1))
        slitlet2d_dn = restoration.denoise_nl_means(slitlet2d_rs,
                                                    patch_size=3,
                                                    patch_distance=2,
                                                    multichannel=False)
        slitlet2d_dn = rescale_array_from_z1z2(slitlet2d_dn, coef_rs)

        # compute basic statistics
        q25, q50, q75 = np.percentile(slitlet2d_dn, q=[25.0, 50.0, 75.0])
        sigmag = 0.7413 * (q75 - q25)  # robust standard deviation
        if abs(self.debugplot) >= 10:
            q16, q84 = np.percentile(slitlet2d_dn, q=[15.87, 84.13])
            print('>>> q16...:', q16)
            print('>>> q25...:', q25)
            print('>>> q50...:', q50)
            print('>>> q75...:', q75)
            print('>>> q84...:', q84)
            print('>>> sigmaG:', sigmag)
        if abs(self.debugplot) in [21, 22]:
            # display initial image with zscale cuts
            title = "Slitlet#" + str(self.islitlet) + \
                    " (locate_unknown_arc_lines, step #1)"
            ximshow(slitlet2d, title=title,
                    first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                    debugplot=self.debugplot)
            # display denoised image with zscale cuts
            title = "Slitlet#" + str(self.islitlet) + \
                    " (locate_unknown_arc_lines, step #2)"
            ximshow(slitlet2d_dn, title=title,
                    first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                    debugplot=self.debugplot)
            # display image with different cuts
            z1z2 = (q50 + times_sigma_threshold * sigmag,
                    q50 + 2 * times_sigma_threshold * sigmag)
            title = "Slitlet#" + str(self.islitlet) + \
                    " (locate_unknown_arc_lines, step #3)"
            ximshow(slitlet2d_dn, title=title, z1z2=z1z2,
                    first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                    debugplot=self.debugplot)

        # determine threshold (using the maximum of q50 + t *sigmag or
        # minimum_threshold)
        threshold = q50 + times_sigma_threshold * sigmag
        if minimum_threshold is not None:
            if minimum_threshold > threshold:
                threshold = minimum_threshold

        # identify objects in slitlet2d above threshold
        labels2d_objects, no_objects = ndimage.label(slitlet2d_dn > threshold)
        if abs(self.debugplot) >= 10:
            print("Number of objects initially found:", no_objects)
        if abs(self.debugplot) in [21, 22]:
            # display all objects identified in the image
            title = "Slitlet#" + str(self.islitlet) + \
                    " (locate_unknown_arc_lines, step #4)"
            z1z2 = (labels2d_objects.min(), labels2d_objects.max())
            ximshow(labels2d_objects, title=title,
                    first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                    cbar_label="Object number",
                    z1z2=z1z2, cmap="nipy_spectral",
                    debugplot=self.debugplot)

        # select arc lines by imposing the criteria based on the
        # dimensions of the detected objects and the intersection with
        # the middle spectrum trail
        slices_possible_arc_lines = ndimage.find_objects(labels2d_objects)
        slices_ok = np.repeat([False], no_objects)  # flag
        for i in range(no_objects):
            if abs(self.debugplot) >= 10:
                print('object', i + 1,
                      '[in np.array coordinates]:',
                      slices_possible_arc_lines[i])
            slice_x = slices_possible_arc_lines[i][1]
            slice_y = slices_possible_arc_lines[i][0]
            # note that the width computation doesn't require to
            # add +1 since slice_x.stop (and slice_y.stop) is
            # already the upper limit +1 (in np.array coordinates)
            delta_x = slice_x.stop - slice_x.start
            delta_y = slice_y.stop - slice_y.start
            # dimensions criterion
            if delta_x <= delta_x_max and delta_y >= delta_y_min:
                # intersection with middle spectrum trail criterion;
                # note that slice_x and slice_y are given in np.array
                # coordinates and are transformed into image coordinates;
                # in addition, -0.5 shift the origin to the lower left
                # corner of the pixel
                xini_slice = slice_x.start + self.bb_nc1_orig - 0.5
                xmiddle_slice = xini_slice + delta_x / 2
                polydum = \
                    self.list_spectrails[self.i_middle_spectrail].poly_funct
                ymiddle_slice = polydum(xmiddle_slice)
                yini_slice = slice_y.start + self.bb_ns1_orig - 0.5
                yend_slice = yini_slice + delta_y
                if yini_slice + min_dist_from_middle <= ymiddle_slice <= \
                        yend_slice - min_dist_from_middle:
                    slices_ok[i] = True

        # generate list with ID of arc lines (note that first object is
        # number 0 and not 1)
        list_slices_ok = []
        for i in range(no_objects):
            if slices_ok[i]:
                list_slices_ok.append(i + 1)
        number_arc_lines = len(list_slices_ok)
        if abs(self.debugplot) >= 10:
            print("\nNumber of arc lines initially identified is:",
                  number_arc_lines)
            if number_arc_lines > 0:
                print("Slice ID of lines passing the selection:\n",
                      list_slices_ok)
        if number_arc_lines == 0:
            return

        # display arc lines
        if abs(self.debugplot) in [21, 22]:
            # display all objects identified in the image
            title = "Slitlet#" + str(self.islitlet) + \
                    " (locate_unknown_arc_lines, step #5)"
            z1z2 = (labels2d_objects.min(),
                    labels2d_objects.max())
            ax = ximshow(labels2d_objects, show=False, title=title,
                         first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                         cbar_label="Object number",
                         z1z2=z1z2, cmap="nipy_spectral",
                         debugplot=self.debugplot)
            # plot rectangle around identified arc lines
            for i in range(no_objects):
                if slices_ok[i]:
                    slice_x = slices_possible_arc_lines[i][1]
                    slice_y = slices_possible_arc_lines[i][0]
                    # note that slice_x and slice_y are given in np.array
                    # coordinates and are transformed into image coordinates;
                    # in addition, -0.5 shift the origin to the lower left
                    # corner of the pixel
                    xini_slice = slice_x.start + self.bb_nc1_orig - 0.5
                    yini_slice = slice_y.start + self.bb_ns1_orig - 0.5
                    # note that the width computation doesn't require to
                    # add +1 since slice_x.stop (and slice_y.stop) is
                    # already the upper limit +1 (in np.array coordinates)
                    xwidth_slice = slice_x.stop - slice_x.start
                    ywidth_slice = slice_y.stop - slice_y.start
                    rect = Rectangle((xini_slice, yini_slice),
                                     xwidth_slice, ywidth_slice,
                                     edgecolor='w', facecolor='none')
                    ax.add_patch(rect)
            # show plot
            pause_debugplot(self.debugplot, pltshow=True)

        # adjust individual arc lines passing the initial selection
        self.list_arc_lines = []  # list of ArcLines
        for k in range(number_arc_lines):  # fit each arc line
            # select points to be fitted for a particular arc line
            xy_tmp = np.where(labels2d_objects == list_slices_ok[k])
            x_tmp = xy_tmp[1] + self.bb_nc1_orig  # use image coordinates
            y_tmp = xy_tmp[0] + self.bb_ns1_orig  # use image coordinates
            w_tmp = slitlet2d_dn[xy_tmp]
            # declare new ArcLine instance
            arc_line = ArcLine()
            # define new ArcLine using a weighted fit
            # (note that it must be X vs Y)
            arc_line.fit(x=x_tmp, y=y_tmp, deg=1, w=w_tmp, y_vs_x=False)
            if len(arc_line.poly_funct.coef) == 2:
                # update list with identified ArcLines
                self.list_arc_lines.append(arc_line)
            else:
                # ignore (sometimes the arc_line.fit returns a constant!)
                pass

        # recompute number_arc_lines just in case in the previous fits
        # some lines have just given a zero degree polynomial
        number_arc_lines = len(self.list_arc_lines)

        # remove arc lines with unexpected slopes
        yfit = np.array([self.list_arc_lines[k].poly_funct.coef[1]
                         for k in range(number_arc_lines)])
        xfit = np.zeros(number_arc_lines)
        # intersection between middle spectrum trail and arc line
        for k in range(number_arc_lines):
            arcline = self.list_arc_lines[k]
            xfit[k], ydum = intersection_spectrail_arcline(
                self.list_spectrails[self.i_middle_spectrail], arcline
            )

        # fit slope versus x-coordinate of the intersection of the arc line
        # with the middle spectrum trail
        if len(yfit) > 5:
            degeff = 5
        else:
            degeff = len(yfit) - 1
        polydum, residum, rejected = polfit_residuals_with_sigma_rejection(
            x=xfit, y=yfit, deg=degeff, times_sigma_reject=4.0,
            xlabel='arc line center (islitlet #' + str(self.islitlet) + ')',
            ylabel='arc line slope', debugplot=0
        )
        # remove rejected arc lines
        if len(rejected) > 0:
            if abs(self.debugplot) >= 10:
                print('Rejecting', sum(rejected),
                      'arc lines with suspicious slopes: Slice ID',
                      [list_slices_ok[k] for k in range(number_arc_lines)
                       if rejected[k]])
            self.list_arc_lines = \
                [self.list_arc_lines[k] for k in range(number_arc_lines)
                 if not rejected[k]]
            # recompute number of arc lines
            number_arc_lines = len(self.list_arc_lines)
            if abs(self.debugplot) >= 10:
                print("\nNumber of arc lines finally identified is:",
                      number_arc_lines)

        if abs(self.debugplot) >= 20:
            # print list of arc lines
            print('\nlist_arc_lines:')
            for k in range(number_arc_lines):
                print(k, '->', self.list_arc_lines[k], '\n')

        # display results
        if abs(self.debugplot) in [21, 22]:
            # generate mask with all the arc-line points passing the selection
            mask_arc_lines = np.zeros_like(slitlet2d_dn)
            for k in list_slices_ok:
                mask_arc_lines[labels2d_objects == k] = 1
            # compute image with only the arc lines passing the selection
            labels2d_arc_lines = labels2d_objects * mask_arc_lines
            # display background image with filtered arc lines
            title = "Slitlet#" + str(self.islitlet) + \
                    " (locate_unknown_arc_lines, step #6)"
            z1z2 = (labels2d_arc_lines.min(),
                    labels2d_arc_lines.max())
            ax = ximshow(labels2d_arc_lines, show=False,
                         first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                         cbar_label="Object number",
                         title=title, z1z2=z1z2, cmap="nipy_spectral",
                         debugplot=self.debugplot)
            # plot weighted fit for each arc line (note that the fit is
            # X vs Y)
            for k in range(number_arc_lines):
                xpol, ypol = self.list_arc_lines[k].linspace_pix()
                ax.plot(xpol, ypol, 'g--')
            # display lower and upper points of each arc line
            x_tmp = [arc_line.xlower_line for arc_line in self.list_arc_lines]
            y_tmp = [arc_line.ylower_line for arc_line in self.list_arc_lines]
            ax.plot(x_tmp, y_tmp, 'w+')
            x_tmp = [arc_line.xupper_line for arc_line in self.list_arc_lines]
            y_tmp = [arc_line.yupper_line for arc_line in self.list_arc_lines]
            ax.plot(x_tmp, y_tmp, 'w+')
            # show plot
            pause_debugplot(self.debugplot, pltshow=True)

    def xy_spectrail_arc_intersections(self, slitlet2d=None):
        """Compute intersection points of spectrum trails with arc lines.

        The member list_arc_lines is updated with new keyword:keyval
        values for each arc line.

        Parameters
        ----------
        slitlet2d : numpy array
            Slitlet image to be displayed with the computed boundaries
            and intersecting points overplotted. This argument is
            optional.

        """

        # protections
        if self.list_arc_lines is None:
            raise ValueError("Arc lines not sought")
        number_spectrum_trails = len(self.list_spectrails)
        if number_spectrum_trails == 0:
            raise ValueError("Number of available spectrum trails is 0")
        number_arc_lines = len(self.list_arc_lines)
        if number_arc_lines == 0:
            raise ValueError("Number of available arc lines is 0")

        # intersection of the arc lines with the spectrum trails
        # (note: the coordinates are computed using pixel values,
        #  ranging from 1 to EMIR_NAXIS1, as given in the original
        #  image reference system ---not in the slitlet image reference
        #  system---)
        self.x_inter_rect = np.array([])  # rectified image coordinates
        self.y_inter_rect = np.array([])  # rectified image coordinates
        for arcline in self.list_arc_lines:
            # middle spectrum trail
            spectrail = self.list_spectrails[self.i_middle_spectrail]
            xroot, yroot = intersection_spectrail_arcline(
                spectrail=spectrail, arcline=arcline
            )
            arcline.x_rectified = xroot
            self.x_inter_rect = np.append(
                self.x_inter_rect, [xroot] * number_spectrum_trails
            )
            for spectrail in self.list_spectrails:
                self.y_inter_rect = np.append(
                    self.y_inter_rect, spectrail.y_rectified)
        #
        self.x_inter_orig = np.array([])  # original image coordinates
        self.y_inter_orig = np.array([])  # original image coordinates
        for arcline in self.list_arc_lines:
            for spectrail in self.list_spectrails:
                xroot, yroot = intersection_spectrail_arcline(
                    spectrail=spectrail, arcline=arcline
                )
                self.x_inter_orig = np.append(self.x_inter_orig, xroot)
                self.y_inter_orig = np.append(self.y_inter_orig, yroot)

        # display intersection points
        if abs(self.debugplot % 10) != 0 and slitlet2d is not None:
            # display image with zscale cuts
            title = "Slitlet#" + str(self.islitlet) + \
                    " (xy_spectrail_arc_intersections)"
            ax = ximshow(slitlet2d, title=title,
                         first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                         show=False)
            # spectrum trails
            for spectrail in self.list_spectrails:
                xdum, ydum = spectrail.linspace_pix(start=self.bb_nc1_orig,
                                                    stop=self.bb_nc2_orig)
                ax.plot(xdum, ydum, 'g')
            # arc lines
            for arcline in self.list_arc_lines:
                xdum, ydum = arcline.linspace_pix(start=self.bb_ns1_orig,
                                                  stop=self.bb_ns2_orig)
                ax.plot(xdum, ydum, 'g')
            # intersection points
            ax.plot(self.x_inter_orig, self.y_inter_orig, 'co')
            ax.plot(self.x_inter_rect, self.y_inter_rect, 'bo')
            # show plot
            pause_debugplot(self.debugplot, pltshow=True)

    def estimate_tt_to_rectify(self, order, slitlet2d=None):
        """Estimate the polynomial transformation to rectify the image.

        Parameters
        ----------
        order : int
            Order of the polynomial transformation.
        slitlet2d : numpy array
            Slitlet image to be displayed with the computed boundaries
            and intersecting points overplotted. This argument is
            optional.

        """

        # protections
        if self.x_inter_orig is None \
                or self.y_inter_orig is None \
                or self.x_inter_rect is None \
                or self.y_inter_rect is None:
            raise ValueError('Intersection points not computed')

        npoints = len(self.x_inter_orig)
        if len(self.y_inter_orig) != npoints \
                or len(self.x_inter_rect) != npoints \
                or len(self.y_inter_rect) != npoints:
            raise ValueError('Unexpected different number of points')

        # IMPORTANT: correct coordinates from origin in order to manipulate
        # coordinates corresponding to image indices
        x_inter_orig_shifted = self.x_inter_orig - self.bb_nc1_orig
        y_inter_orig_shifted = self.y_inter_orig - self.bb_ns1_orig
        x_inter_rect_shifted = self.x_inter_rect - self.bb_nc1_orig
        y_inter_rect_shifted = self.y_inter_rect - self.bb_ns1_orig

        # compute 2D transformation
        self.ttd_order = order
        self.ttd_aij, self.ttd_bij = compute_distortion(
            x_inter_orig_shifted, y_inter_orig_shifted,
            x_inter_rect_shifted, y_inter_rect_shifted,
            order,
            self.debugplot
        )
        self.tti_aij, self.tti_bij = compute_distortion(
            x_inter_rect_shifted, y_inter_rect_shifted,
            x_inter_orig_shifted, y_inter_orig_shifted,
            order,
            self.debugplot
        )

        # display slitlet with intersection points and grid indicating
        # the fitted transformation
        if abs(self.debugplot % 10) != 0 and slitlet2d is not None:
            # display image with zscale cuts
            title = "Slitlet#" + str(self.islitlet) + \
                    " (estimate_tt_to_rectify)"
            ax = ximshow(slitlet2d, title=title,
                         first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                         show=False)
            # intersection points
            ax.plot(self.x_inter_orig, self.y_inter_orig, 'co')
            ax.plot(self.x_inter_rect, self.y_inter_rect, 'bo')
            # grid with fitted transformation: spectrum trails
            xx = np.arange(0, self.bb_nc2_orig - self.bb_nc1_orig + 1,
                           dtype=np.float)
            for spectrail in self.list_spectrails:
                yy0 = spectrail.y_rectified
                yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
                ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b")
                xxx, yyy = fmap(self.ttd_order, self.ttd_aij, self.ttd_bij,
                                xx, yy)
                ax.plot(xxx + self.bb_nc1_orig, yyy + self.bb_ns1_orig, "g")
            # grid with fitted transformation: arc lines
            ylower_line = \
                self.list_spectrails[self.i_lower_spectrail].y_rectified
            yupper_line = \
                self.list_spectrails[self.i_upper_spectrail].y_rectified
            n_points = int(yupper_line - ylower_line + 0.5) + 1
            yy = np.linspace(ylower_line - self.bb_ns1_orig,
                             yupper_line - self.bb_ns1_orig,
                             num=n_points,
                             dtype=np.float)
            for arc_line in self.list_arc_lines:
                xline = arc_line.x_rectified - self.bb_nc1_orig
                xx = np.array([xline] * n_points)
                ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b")
                xxx, yyy = fmap(self.ttd_order, self.ttd_aij, self.ttd_bij,
                                xx, yy)
                ax.plot(xxx + self.bb_nc1_orig, yyy + self.bb_ns1_orig, "c")
            # show plot
            pause_debugplot(self.debugplot, pltshow=True)

    def rectify(self, slitlet2d, resampling, transformation, inverse=False):
        """Rectify slitlet using computed transformation.

        Parameters
        ----------
        slitlet2d : numpy array
            Image containing the 2d slitlet image.
        resampling : int
            1: nearest neighbour, 2: flux preserving interpolation.
        transformation : int
            1: initial, 2: modeled
        inverse : bool
            If true, the inverse rectification transformation is
            employed.

        Returns
        -------
        slitlet2d_rect : numpy array
            Rectified slitlet image.

        """

        if resampling not in [1, 2]:
            raise ValueError("Unexpected resampling value=" + str(resampling))

        if transformation not in [1, 2]:
            raise ValueError("Unexpected transformation value=" +
                             str(transformation))

        # verify image dimension
        naxis2, naxis1 = slitlet2d.shape
        if naxis1 != self.bb_nc2_orig - self.bb_nc1_orig + 1:
            raise ValueError("Unexpected slitlet2d_rect naxis1")
        if naxis2 != self.bb_ns2_orig - self.bb_ns1_orig + 1:
            raise ValueError("Unexpected slitlet2d_rect naxis2")

        # transformation to be employed (and direction)
        if transformation == 1:
            if inverse:
                aij = self.tti_aij
                bij = self.tti_bij
            else:
                aij = self.ttd_aij
                bij = self.ttd_bij
        else:
            if inverse:
                aij = self.tti_aij_longslit_model
                bij = self.tti_bij_longslit_model
            else:
                aij = self.ttd_aij_longslit_model
                bij = self.ttd_bij_longslit_model

        # rectify image
        slitlet2d_rect = rectify2d(
            image2d=slitlet2d,
            aij=aij,
            bij=bij,
            resampling=resampling
        )

        if abs(self.debugplot % 10) != 0:
            title = "Slitlet#" + str(self.islitlet) + " (rectify)"
            ax = ximshow(slitlet2d_rect, title=title,
                         first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                         show=False)
            if self.list_arc_lines is not None:
                # intersection points
                ax.plot(self.x_inter_rect, self.y_inter_rect, 'bo')
            # grid with fitted transformation: spectrum trails
            xx = np.arange(0, self.bb_nc2_orig - self.bb_nc1_orig + 1,
                           dtype=np.float)
            for spectrail in self.list_spectrails:
                yy0 = spectrail.y_rectified
                yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
                ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b")
            for spectrail in self.list_frontiers:
                yy0 = spectrail.y_rectified
                yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
                ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b:")
            # grid with fitted transformation: arc lines
            ylower_line = \
                self.list_spectrails[self.i_lower_spectrail].y_rectified
            yupper_line = \
                self.list_spectrails[self.i_upper_spectrail].y_rectified
            n_points = int(yupper_line - ylower_line + 0.5) + 1
            yy = np.linspace(ylower_line - self.bb_ns1_orig,
                             yupper_line - self.bb_ns1_orig,
                             num=n_points,
                             dtype=np.float)
            if self.list_arc_lines is not None:
                for arc_line in self.list_arc_lines:
                    xline = arc_line.x_rectified - self.bb_nc1_orig
                    xx = np.array([xline] * n_points)
                    ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b")
            # show plot
            pause_debugplot(self.debugplot, pltshow=True)

        return slitlet2d_rect

    def median_spectrum_from_rectified_image(self, slitlet2d_rect,
                                             sigma_gaussian_filtering=0,
                                             nwinwidth_initial=5,
                                             nwinwidth_refined=7,
                                             times_sigma_threshold=5,
                                             minimum_threshold=None,
                                             npix_avoid_border=0,
                                             nbrightlines=None):
        """Median spectrum and line peaks from rectified image.

        In order to avoid the line ghosts, the line peaks are identified
        independently in the upper and lower halves of the rectified
        image. The final peaks correspond to lines that appear in both
        spectra.

        Parameters
        ----------
        slitlet2d_rect : numpy array
            Rectified slitlet image.
        sigma_gaussian_filtering : float
            Sigma of the gaussian filter to be applied to the spectrum
            in order to avoid problems with saturated lines. This
            filtering is skipped when this parameter is <= 0.
        nwinwidth_initial : int
            Width of the window where each peak must be found using
            the initial method (approximate)
        nwinwidth_refined : int
            Width of the window where each peak location will be
            refined.
        times_sigma_threshold : float
            Times (robust) sigma above the median of the image to set
            the minimum threshold when searching for line peaks.
        minimum_threshold : float or None
            Minimum value of the threshold.
        npix_avoid_border : int
            Number of pixels at the borders of the spectrum where peaks
            are not considered. If zero, the actual number will be
            given by nwinwidth_initial.
        nbrightlines : list
            List with maximum number of brightest lines to be employed
            in the wavelength calibration. The length of the list
            indicates the number of equal-size subintervals in the
            wavelength calibration direction to be considered.
            If this value is [0] or None, all the detected lines will be
            employed.

        Returns
        -------
        sp0 : numpy array
            Median spectrum.
        fxpeaks : numpy array
            Refined location of arc lines (in array index scale).

        """

        # protections
        naxis2, naxis1 = slitlet2d_rect.shape
        if naxis1 != self.bb_nc2_orig - self.bb_nc1_orig + 1:
            raise ValueError("Unexpected slitlet2d_rect naxis1")
        if naxis2 != self.bb_ns2_orig - self.bb_ns1_orig + 1:
            raise ValueError("Unexpected slitlet2d_rect naxis2")

        # lower, middle and upper spectrum trails
        ylower_line = \
            self.list_spectrails[self.i_lower_spectrail].y_rectified
        ymiddle_line = \
            self.list_spectrails[self.i_middle_spectrail].y_rectified
        yupper_line = \
            self.list_spectrails[self.i_upper_spectrail].y_rectified

        ilower = int(ylower_line + 0.5) - self.bb_ns1_orig
        imiddle = int(ymiddle_line + 0.5) - self.bb_ns1_orig
        iupper = int(yupper_line + 0.5) - self.bb_ns1_orig

        # median spectra using different image regions
        sp0_ini = np.median(slitlet2d_rect[ilower:(iupper + 1), :], axis=0)
        sp1_ini = np.median(slitlet2d_rect[ilower:(imiddle + 1), :], axis=0)
        sp2_ini = np.median(slitlet2d_rect[imiddle:(iupper + 1), :], axis=0)

        # gaussian filtering when requested (to avoid line saturation)
        if sigma_gaussian_filtering > 0:
            sp0 = ndimage.filters.gaussian_filter(
                sp0_ini,
                sigma=sigma_gaussian_filtering
            )
            sp1 = ndimage.filters.gaussian_filter(
                sp1_ini,
                sigma=sigma_gaussian_filtering
            )
            sp2 = ndimage.filters.gaussian_filter(
                sp2_ini,
                sigma=sigma_gaussian_filtering
            )
        else:
            sp0 = np.copy(sp0_ini)
            sp1 = np.copy(sp1_ini)
            sp2 = np.copy(sp2_ini)

        # compute threshold
        q25, q50, q75 = np.percentile(sp0, q=[25.0, 50.0, 75.0])
        sigma_g = 0.7413 * (q75 - q25)  # robust standard deviation
        threshold = q50 + times_sigma_threshold * sigma_g
        if abs(self.debugplot) >= 10:
            print("median...........:", q50)
            print("robuts std.......:", sigma_g)
            print("threshold........:", threshold)
        if minimum_threshold is not None:
            if minimum_threshold > threshold:
                threshold = minimum_threshold
        if abs(self.debugplot) >= 10:
            print("minimum threshold:", minimum_threshold)
            print("final threshold..:", threshold)

        # initial location of the peaks (integer values)
        ixpeaks0 = find_peaks_spectrum(sp0, nwinwidth=nwinwidth_initial,
                                       threshold=threshold,
                                       debugplot=self.debugplot)

        # peaks in the lower and upper regions
        ixpeaks1 = find_peaks_spectrum(sp1, nwinwidth=nwinwidth_initial,
                                       threshold=threshold,
                                       debugplot=self.debugplot)
        ixpeaks2 = find_peaks_spectrum(sp2, nwinwidth=nwinwidth_initial,
                                       threshold=threshold,
                                       debugplot=self.debugplot)

        # the peaks are valid if the are also found in the lower and
        # upper regions (with a tolerance of +1 or -1 pixel)
        ixpeaks = []
        for ixpeak in ixpeaks0:
            l1 = ixpeak in np.concatenate((ixpeaks1, ixpeaks1+1, ixpeaks1-1))
            l2 = ixpeak in np.concatenate((ixpeaks2, ixpeaks2+1, ixpeaks2-1))
            if l1 and l2:
                ixpeaks.append(ixpeak)
        ixpeaks = np.array(ixpeaks)
        if abs(self.debugplot) >= 10:
            print("Merged initial list of peaks:\n", ixpeaks)

        # remove peaks too close to any of the borders of the spectrum
        if npix_avoid_border > 0:
            lok_ini = ixpeaks >= npix_avoid_border
            lok_end = ixpeaks <= len(sp0) - 1 - npix_avoid_border
            ixpeaks = ixpeaks[lok_ini * lok_end]

        # select a maximum number of brightest lines in each region
        if nbrightlines is None:
            pass
        elif len(nbrightlines) == 1 and nbrightlines[0] == 0:
            pass
        else:
            if abs(self.debugplot) >= 10:
                print('nbrightlines =', nbrightlines)
                print('ixpeaks in whole spectrum:\n', ixpeaks)
            region_size = (naxis1-1)/len(nbrightlines)
            ixpeaks_filtered = np.array([], dtype=int)
            for iregion, nlines_in_region in enumerate(nbrightlines):
                if nlines_in_region > 0:
                    imin = int(iregion * region_size)
                    imax = int((iregion + 1) * region_size)
                    if iregion > 0:
                        imin += 1
                    ixpeaks_region = \
                        ixpeaks[np.logical_and(ixpeaks >= imin,
                                               ixpeaks <= imax)]
                    if len(ixpeaks_region) > 0:
                        peak_fluxes = sp0[ixpeaks_region]
                        spos = peak_fluxes.argsort()
                        ixpeaks_tmp = ixpeaks_region[spos[-nlines_in_region:]]
                        ixpeaks_tmp.sort()  # in-place sort
                        if abs(self.debugplot) >= 10:
                            print('ixpeaks in region........:\n', ixpeaks_tmp)
                        ixpeaks_filtered = np.concatenate((ixpeaks_filtered,
                                                           ixpeaks_tmp))
            ixpeaks = ixpeaks_filtered
            if abs(self.debugplot) >= 10:
                print('ixpeaks filtered.........:\n', ixpeaks)

        # refined location of the peaks (float values)
        fxpeaks, sxpeaks = refine_peaks_spectrum(sp0, ixpeaks,
                                                 nwinwidth=nwinwidth_refined,
                                                 method="gaussian")

        if abs(self.debugplot) % 10 != 0:
            x = np.arange(self.bb_nc1_orig, self.bb_nc2_orig + 1)
            title = "Slitlet#" + str(self.islitlet) + " (median spectrum)"
            ax = ximplotxy(x, sp1, show=False, title=title,
                           xlabel='pixel coordinate (from 1 to NAXIS1)',
                           ylabel='number of counts',
                           **{'marker': ' ', 'label': 'lower region'})
            ax.plot(x, sp2, label='upper region')
            ax.plot(x, sp0, label='whole region')
            # mark peak location
            ax.plot(ixpeaks + self.bb_nc1_orig,
                    sp0[ixpeaks], 'o', label="initial location")
            ax.plot(fxpeaks + self.bb_nc1_orig,
                    sp0[ixpeaks], 'o', label="refined location")
            ax.legend()
            pause_debugplot(self.debugplot, pltshow=True, tight_layout=False)

        # return median spectrum and refined peak location
        return sp0, fxpeaks


class Slitlet2D(object):
    """Slitlet2D class definition.

    Parameters
    ----------
    islitlet : int
        Slitlet number.
    rectwv_coeff : instance of RectWaveCoeff
        Object storing the JSON input file where the
        the rectification and wavelength calibration transformations
        for a particular instrument configuration are stored.
    debugplot : int
        Debugging level for messages and plots. For details see
        'numina.array.display.pause_debugplot.py'.

    Attributes
    ----------
    islitlet : int
        Slitlet number.
    csu_bar_left : float
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
    bb_nc1_orig : int
        Minimum X coordinate of the enclosing bounding box (in pixel
        units) in the original image.
    bb_nc2_orig : int
        Maximum X coordinate of the enclosing bounding box (in pixel
        units) in the original image.
    bb_ns1_orig : int
        Minimum Y coordinate of the enclosing bounding box (in pixel
        units) in the original image.
    bb_ns2_orig : int
        Maximum Y coordinate of the enclosing bounding box (in pixel
        units) in the original image.
    x0_reference : float
        X coordinate where the rectified y0_reference_middle is computed
        as the Y coordinate of the middle spectrum trail. The same value
        is used for all the available spectrum trails.
    list_spectrails: list of numpy.polynomial.Polynomial instances
        List of spectrum trails defined (lower, middle and upper).
    list_frontiers: list of numpy.polynomial.Polynomial instances
        List of spectrum trails defining the slitlet frontiers (lower
        and upper).
    y0_reference_lower: float
        Y coordinate corresponding to the lower spectrum trail computed
        at x0_reference. This value is employed as the Y coordinate of
        the lower spectrum trail of the rectified slitlet.
    y0_reference_middle: float
        Y coordinate corresponding to the middle spectrum trail computed
        at x0_reference. This value is employed as the Y coordinate of
        the middle spectrum trail of the rectified slitlet.
    y0_reference_upper: float
        Y coordinate corresponding to the upper spectrum trail computed
        at x0_reference. This value is employed as the Y coordinate of
        the upper spectrum trail of the rectified slitlet.
    y0_frontier_lower: float
        Y coordinate corresponding to the lower frontier computed at
        x0_reference.
    y0_frontier_upper: float
        Y coordinate corresponding to the upper frontier computed at
        x0_reference.
    ttd_order : int or None
        Polynomial order corresponding to the rectification
        transformation.
    ttd_aij : numpy array
        Polynomial coefficents corresponding to the direct
        rectification transformation coefficients a_ij.
    ttd_bij : numpy array
        Polynomial coefficents corresponding to the direct
        rectification transformation coefficients b_ij.
    tti_aij : numpy array
        Polynomial coefficents corresponding to the inverse
        rectification transformation coefficients a_ij.
    tti_bij : numpy array
        Polynomial coefficents corresponding to the inverse
        rectification transformation coefficients b_ij.
    wpoly : Polynomial instance
        Wavelength calibration polynomial, providing the
        wavelength as a function of pixel number (running from 1 to
        NAXIS1).
    debugplot : int
        Debugging level for messages and plots. For details see
        'numina.array.display.pause_debugplot.py'.


    """

    def __init__(self, islitlet, rectwv_coeff, debugplot):
        # slitlet number
        self.islitlet = islitlet

        tmpcontent = rectwv_coeff.contents[islitlet-1]

        # csu configuration
        self.csu_bar_left = tmpcontent['csu_bar_left']
        self.csu_bar_right = tmpcontent['csu_bar_right']
        self.csu_bar_slit_center = tmpcontent['csu_bar_slit_center']
        self.csu_bar_slit_width = tmpcontent['csu_bar_slit_width']

        # horizontal and vertical bounding box
        self.bb_nc1_orig = tmpcontent['bb_nc1_orig']
        self.bb_nc2_orig = tmpcontent['bb_nc2_orig']
        self.bb_ns1_orig = tmpcontent['bb_ns1_orig']
        self.bb_ns2_orig = tmpcontent['bb_ns2_orig']

        # reference abscissa
        self.x0_reference = tmpcontent['x0_reference']

        # list of spectrum trails (lower, middle, and upper)
        self.list_spectrails = []
        for idum, cdum in zip(range(3), ['lower', 'middle', 'upper']):
            coeff = tmpcontent['spectrail']['poly_coef_' + cdum]
            self.list_spectrails.append(np.polynomial.Polynomial(coeff))

        # define reference ordinates using lower, middle and upper spectrails
        # evaluated at x0_reference
        self.y0_reference_lower = tmpcontent['y0_reference_lower']
        self.y0_reference_middle = tmpcontent['y0_reference_middle']
        self.y0_reference_upper = tmpcontent['y0_reference_upper']

        # list of frontiers (lower and upper)
        self.list_frontiers = []
        for idum, cdum in zip(range(2), ['lower', 'upper']):
            coeff = tmpcontent['frontier']['poly_coef_' + cdum]
            self.list_frontiers.append(np.polynomial.Polynomial(coeff))

        # define frontier ordinates at x0_reference
        self.y0_frontier_lower = tmpcontent['y0_frontier_lower']
        self.y0_frontier_upper = tmpcontent['y0_frontier_upper']

        # Rectification coefficients
        self.ttd_aij = tmpcontent['ttd_aij']
        self.ttd_bij = tmpcontent['ttd_bij']
        self.tti_aij = tmpcontent['tti_aij']
        self.tti_bij = tmpcontent['tti_bij']
        # determine order from number of coefficients
        ncoef = len(self.ttd_aij)
        self.ttd_order = order_fmap(ncoef)

        # Wavelength calibration coefficients
        self.wpoly = tmpcontent['wpoly_coeff']

        # debugplot
        self.debugplot = debugplot

    def __repr__(self):
        """Define printable representation of a Slitlet2D instance."""

        # string with all the information
        output = "<Slilet2D instance>\n" + \
            "- islitlet...........: " + \
                 str(self.islitlet) + "\n" + \
            "- csu_bar_left.......: " + \
                 str(self.csu_bar_left) + "\n" + \
            "- csu_bar_right......: " + \
                 str(self.csu_bar_right) + "\n" + \
            "- csu_bar_slit_center: " + \
                 str(self.csu_bar_slit_center) + "\n" + \
            "- csu_bar_slit_width.: " + \
                 str(self.csu_bar_slit_width) + "\n" + \
            "- x0_reference.......: " + \
                 str(self.x0_reference) + "\n" + \
            "- y0_reference_lower.: " + \
                 str(self.y0_reference_lower) + "\n" + \
            "- y0_reference_middle: " + \
                 str(self.y0_reference_middle) + "\n" + \
            "- y0_reference_upper.: " + \
                 str(self.y0_reference_upper) + "\n" + \
            "- y0_frontier_lower..: " + \
                 str(self.y0_frontier_lower) + "\n" + \
            "- y0_frontier_upper..: " + \
                 str(self.y0_frontier_upper) + "\n" + \
            "- bb_nc1_orig........: " + \
                 str(self.bb_nc1_orig) + "\n" + \
            "- bb_nc2_orig........: " + \
                 str(self.bb_nc2_orig) + "\n" + \
            "- bb_ns1_orig........: " + \
                 str(self.bb_ns1_orig) + "\n" + \
            "- bb_ns2_orig........: " + \
                 str(self.bb_ns2_orig) + "\n" + \
            "- lower spectrail....:\n\t" + \
                 str(self.list_spectrails[0]) + "\n" + \
            "- middle spectrail...:\n\t" + \
                 str(self.list_spectrails[1]) + "\n" + \
            "- upper spectrail....:\n\t" + \
                 str(self.list_spectrails[2]) + "\n" + \
            "- lower frontier.....:\n\t" + \
                 str(self.list_frontiers[0]) + "\n" + \
            "- upper frontier.....:\n\t" + \
                 str(self.list_frontiers[1]) + "\n" + \
            "- ttd_order..........: " + str(self.ttd_order) + "\n" + \
            "- ttd_aij............:\n\t" + str(self.ttd_aij) + "\n" + \
            "- ttd_bij............:\n\t" + str(self.ttd_bij) + "\n" + \
            "- tti_aij............:\n\t" + str(self.tti_aij) + "\n" + \
            "- tti_bij............:\n\t" + str(self.tti_bij) + "\n" + \
            "- wpoly..............:\n\t" + str(self.wpoly) + "\n" + \
            "- debugplot...................: " + \
            str(self.debugplot)

        return output

    def extract_slitlet2d(self, image_2k2k):
        """Extract slitlet 2d image from image with original EMIR dimensions.

        Parameters
        ----------
        image_2k2k : numpy array
            Original image (dimensions EMIR_NAXIS1 * EMIR_NAXIS2)

        Returns
        -------
        slitlet2d : numpy array
            Image corresponding to the slitlet region defined by its
            bounding box.

        """

        # protections
        naxis2, naxis1 = image_2k2k.shape
        if naxis1 != EMIR_NAXIS1:
            raise ValueError('Unexpected naxis1')
        if naxis2 != EMIR_NAXIS2:
            raise ValueError('Unexpected naxis2')

        # extract slitlet region
        slitlet2d = image_2k2k[(self.bb_ns1_orig - 1):self.bb_ns2_orig,
                               (self.bb_nc1_orig - 1):self.bb_nc2_orig]

        # transform to float
        slitlet2d = slitlet2d.astype(np.float)

        # display slitlet2d with boundaries and middle spectrum trail
        if abs(self.debugplot) in [21, 22]:
            self.ximshow_unrectified(slitlet2d)

        # return slitlet image
        return slitlet2d

    def rectify(self, slitlet2d, resampling, inverse=False):
        """Rectify slitlet using computed transformation.

        Parameters
        ----------
        slitlet2d : numpy array
            Image containing the 2d slitlet image.
        resampling : int
            1: nearest neighbour, 2: flux preserving interpolation.
        inverse : bool
            If true, the inverse rectification transformation is
            employed.

        Returns
        -------
        slitlet2d_rect : numpy array
            Rectified slitlet image.

        """

        if resampling not in [1, 2]:
            raise ValueError("Unexpected resampling value=" + str(resampling))

        # check image dimension
        naxis2, naxis1 = slitlet2d.shape
        if naxis1 != self.bb_nc2_orig - self.bb_nc1_orig + 1:
            raise ValueError("Unexpected slitlet2d_rect naxis1")
        if naxis2 != self.bb_ns2_orig - self.bb_ns1_orig + 1:
            raise ValueError("Unexpected slitlet2d_rect naxis2")

        if inverse:
            aij = self.tti_aij
            bij = self.tti_bij
        else:
            aij = self.ttd_aij
            bij = self.ttd_bij

        # rectify image
        slitlet2d_rect = rectify2d(
            image2d=slitlet2d,
            aij=aij,
            bij=bij,
            resampling=resampling
        )

        if abs(self.debugplot % 10) != 0:
            if inverse:
                self.ximshow_unrectified(slitlet2d_rect)
            else:
                self.ximshow_rectified(slitlet2d_rect)

        return slitlet2d_rect

    def ximshow_unrectified(self, slitlet2d):
        """Display unrectified image with spectrails and frontiers.

        Parameters
        ----------
        slitlet2d : numpy array
            Array containing the unrectified slitlet image.

        """

        title = "Slitlet#" + str(self.islitlet)
        ax = ximshow(slitlet2d, title=title,
                     first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                     show=False)
        xdum = np.linspace(1, EMIR_NAXIS1, num=EMIR_NAXIS1)
        ylower = self.list_spectrails[0](xdum)
        ax.plot(xdum, ylower, 'b-')
        ymiddle = self.list_spectrails[1](xdum)
        ax.plot(xdum, ymiddle, 'b--')
        yupper = self.list_spectrails[2](xdum)
        ax.plot(xdum, yupper, 'b-')
        ylower_frontier = self.list_frontiers[0](xdum)
        ax.plot(xdum, ylower_frontier, 'b:')
        yupper_frontier = self.list_frontiers[1](xdum)
        ax.plot(xdum, yupper_frontier, 'b:')
        pause_debugplot(debugplot=self.debugplot, pltshow=True)

    def ximshow_rectified(self, slitlet2d_rect):
        """Display rectified image with spectrails and frontiers.

        Parameters
        ----------
        slitlet2d_rect : numpy array
            Array containing the rectified slitlet image

        """

        title = "Slitlet#" + str(self.islitlet) + " (rectify)"
        ax = ximshow(slitlet2d_rect, title=title,
                     first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                     show=False)
        # grid with fitted transformation: spectrum trails
        xx = np.arange(0, self.bb_nc2_orig - self.bb_nc1_orig + 1,
                       dtype=np.float)
        for spectrail in self.list_spectrails:
            yy0 = spectrail(self.x0_reference)
            yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
            ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b")
        for spectrail in self.list_frontiers:
            yy0 = spectrail(self.x0_reference)
            yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
            ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b:")
        # show plot
        pause_debugplot(self.debugplot, pltshow=True)


def rescale_array_to_z1z2(array, z1z2=(-1.0, 1.0)):
    """Rescale the values in a numpy array to the [z1,z2] interval.

    The transformation is carried out following the relation
    array_rs = b_flux * array - c_flux
    as explained in Appendix B1 of Cardiel (2009, MNRAS, 396, 680)

    Parameters
    ----------
    array : numpy array
        Numpy array to be rescaled.
    z1z2 : tuple, floats
        Minimum and maximum values in the returned array.

    Returns
    -------
    array_rs : numpy array
        Array with rescaled values.
    coef_rs : tuple, floats
        Coefficients b_flux and c_flux employed in the rescaling
        operation.

    """

    if type(array) is not np.ndarray:
        raise ValueError("array=" + str(array) + " must be a numpy.ndarray")

    array_min = array.min()
    array_max = array.max()

    z1, z2 = z1z2
    delta = array_max - array_min
    b_flux = (z2 - z1) / delta
    c_flux = (z2 * array_min - z1 * array_max) / delta

    array_rs = b_flux * array - c_flux

    return array_rs, (b_flux, c_flux)


def rescale_array_from_z1z2(array_rs, coef_rs=None):
    """Restore the values in a numpy array rescaled to the [z1,z2] interval.

    The transformation is carried out following the relation
    array = (array_rs + c_flux)/b_flux
    as explained in Appendix B1 of Cardiel (2009, MNRAS, 396, 680)

    Parameters
    ----------
    array_rs : numpy array
        Numpy array previously rescaled to the [z1,z2] interval with
        the function rescale_array_to_z1z2().
    coef_rs : tuple, floats
        Coefficients b_flux and c_flux previously employed in the
        rescaling operation. This tuple is one of the parameters
        returned by function_rescale_array_to_z1z2().

    Returns
    -------
    array : numpy array
        Array with restored values.

    """

    if type(array_rs) is not np.ndarray:
        raise ValueError(
            "array_rs=" + str(array_rs) + "must be a numpy.ndarray")

    b_flux, c_flux = coef_rs

    array = (array_rs + c_flux) / b_flux

    return array


def set_wv_parameters(filter_name, grism_name):
    """Set wavelength calibration parameters for rectified images.

    Parameters
    ----------
    filter_name : str
        Filter name.
    grism_name : str
        Grism name.

    Returns
    -------
    wv_parameters : dictionary
        Python dictionary containing relevant wavelength calibration
        parameters:
        - crpix1_enlarged : float
        - crval1_enlarged : float
        - cdelt1_enlarged : float
        - naxis1_enlarged: int
        - islitlet_min : int
          Minimium slitlet number.
        - islitlet_max : int
          Maximium slitlet number.
        - nbrightlines : python list
          List of integers containing the number of brightlines to be
          used in the initial wavelength calibration.
        - poly_crval1_linear : numpy.polynomial.Polynomial instance
          Polynomial providing the value of CRVAL1_linear as a function
          of csu_bar_slit_center.
        - poly_cdelt1_linear : numpy.polynomial.Polynomial instance
          Polynomial providing the value of CDELT1_linear as a function
          of csu_bar_slit_center.

    """

    # protections
    if filter_name not in EMIR_VALID_FILTERS:
        raise ValueError('Unexpected filter_name:', filter_name)
    if grism_name not in EMIR_VALID_GRISMS:
        raise ValueError('Unexpected grism_name:', grism_name)

    # intialize output
    wv_parameters = {}
    # set parameters
    wv_parameters['crpix1_enlarged'] = 1.0
    if grism_name == "J" and filter_name == "J":
        wv_parameters['islitlet_min'] = 2
        wv_parameters['islitlet_max'] = 54
        wv_parameters['nbrightlines'] = [18]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            1.25137094e+04,
            -4.81553731e+00,
            4.70039758e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            7.74133267e-01,
            -4.72423718e-05,
            2.79842624e-08
        ])
        wv_parameters['crval1_enlarged'] = 11200.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 0.77        # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = 3400        # pixels
    elif grism_name == "H" and filter_name == "H":
        wv_parameters['islitlet_min'] = 2
        wv_parameters['islitlet_max'] = 54
        wv_parameters['nbrightlines'] = [15]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            1.65561291e+04,
            -7.62779414e+00,
            7.52228172e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            1.21372302e+00,
            3.77209658e-06,
            -1.07467162e-07
        ])
        wv_parameters['crval1_enlarged']= 14500.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 1.2200      # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = 3400        # pixels
    elif grism_name == "K" and filter_name == "Ksp":
        wv_parameters['islitlet_min'] = 2
        wv_parameters['islitlet_max'] = 54
        wv_parameters['nbrightlines'] = [15]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            2.21095313e+04,
            -1.08900414e+01,
            9.66474839e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            1.72596244e+00,
            2.85573046e-05,
            -1.30027272e-07
        ])
        wv_parameters['crval1_enlarged'] = 19100.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 1.7300      # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = 3400        # pixels
    elif grism_name == "LR" and filter_name == "YJ":
        wv_parameters['islitlet_min'] = 4
        wv_parameters['islitlet_max'] = 55
        wv_parameters['nbrightlines'] = None
        wv_parameters['poly_crval1_linear'] = None
        wv_parameters['poly_cdelt1_linear'] = None
        wv_parameters['crval1_enlarged'] = None        # Angstroms
        wv_parameters['cdelt1_enlarged'] = None        # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = None        # pixels
    elif grism_name == "LR" and filter_name == "HK":
        wv_parameters['islitlet_min'] = 4
        wv_parameters['islitlet_max'] = 55
        wv_parameters['nbrightlines'] = None
        wv_parameters['poly_crval1_linear'] = None
        wv_parameters['poly_cdelt1_linear'] = None
        wv_parameters['crval1_enlarged'] = None        # Angstroms
        wv_parameters['cdelt1_enlarged'] = None        # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = None        # pixels
    else:
        print("filter_name..:", filter_name)
        print("grism_name...:", grism_name)
        raise ValueError("invalid filter_name and grism_name combination")

    return wv_parameters


def islitlet_progress(islitlet, islitlet_max):
    """Auxiliary function to print out progress in loop of slitlets.

    Parameters
    ----------
    islitlet : int
        Current slitlet number.
    islitlet_max : int
        Maximum slitlet number.

    """
    if islitlet % 10 == 0:
        cout = str(islitlet // 10)
    else:
        cout = '.'
    sys.stdout.write(cout)
    if islitlet == islitlet_max:
        sys.stdout.write('\n')
    sys.stdout.flush()


def rectwv_coeff_from_arc_image(reduced_image,
                                bound_param,
                                lines_catalog,
                                args_nbrightlines=None,
                                args_ymargin_bb=2,
                                args_remove_sp_background=True,
                                args_times_sigma_threshold=10,
                                args_order_fmap=2,
                                args_sigma_gaussian_filtering=2,
                                args_margin_npix=50,
                                args_poldeg_initial=3,
                                args_poldeg_refined=5,
                                args_interactive=False,
                                args_threshold_wv=0,
                                args_ylogscale=False,
                                args_pdf=None,
                                args_geometry=(0,0,640,480),
                                debugplot=0):
    """Evaluate rect.+wavecal. coefficients from arc image

    Parameters
    ----------
    reduced_image : HDUList object
        Image with preliminary basic reduction: bpm, bias, dark and
        flatfield.
    bound_param : RefinedBoundaryModelParam instance
        Refined boundary model.
    lines_catalog : Numpy array
        2D numpy array with the contents of the master file with the
        expected arc line wavelengths.
    debugplot : int
            Debugging level for messages and plots. For details see
            'numina.array.display.pause_debugplot.py'.

    Returns
    -------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration of the input arc image.

    """

    logger = logging.getLogger(__name__)

    # header and data array
    header = reduced_image[0].header
    image2d = reduced_image[0].data

    # check grism and filter
    filter_name = header['filter']
    logger.debug('Filter: ' + filter_name)
    if filter_name != bound_param.tags['filter']:
        raise ValueError('Filter name does not match!')
    grism_name = header['grism']
    logger.debug('Grism: ' + grism_name)
    if grism_name != bound_param.tags['grism']:
        raise ValueError('Grism name does not match!')

    # read the CSU configuration from the image header
    csu_conf = CsuConfiguration.define_from_header(header)
    logger.debug(csu_conf)

    # read the DTU configuration from the image header
    dtu_conf = DtuConfiguration.define_from_header(header)
    logger.debug(dtu_conf)

    # set boundary parameters
    parmodel = bound_param.meta_info['parmodel']
    params = bound_params_from_dict(bound_param.__getstate__())
    if abs(debugplot) >= 10:
        print('-' * 83)
        print('* FITTED BOUND PARAMETERS')
        params.pretty_print()
        pause_debugplot(debugplot)

    # determine parameters according to grism+filter combination
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    islitlet_min = wv_parameters['islitlet_min']
    islitlet_max = wv_parameters['islitlet_max']
    if args_nbrightlines is None:
        nbrightlines = wv_parameters['nbrightlines']
    else:
        nbrightlines = [int(idum) for idum in args_nbrightlines.split(',')]
    poly_crval1_linear = wv_parameters['poly_crval1_linear']
    poly_cdelt1_linear = wv_parameters['poly_cdelt1_linear']

    # list of slitlets to be computed
    list_slitlets = list(range(islitlet_min, islitlet_max + 1))
    print('list_slitlets:\n', list_slitlets)

    # read master arc line wavelengths (only brightest lines)
    wv_master = read_wv_master_from_array(
        master_table=lines_catalog, lines='brightest', debugplot=debugplot
    )

    # read master arc line wavelengths (whole data set)
    wv_master_all = read_wv_master_from_array(
        master_table=lines_catalog, lines='all', debugplot=debugplot
    )

    # check that the arc lines in the master file are properly sorted
    # in ascending order
    for i in range(len(wv_master_all) - 1):
        if wv_master_all[i] >= wv_master_all[i + 1]:
            print('>>> wavelengths: ', wv_master_all[i], wv_master_all[i + 1])
            raise ValueError('Arc lines are not sorted in master file')

    # ---

    image2d_55sp = np.zeros((EMIR_NBARS, EMIR_NAXIS1))

    # compute rectification transformation and wavelength calibration
    # polynomials

    measured_slitlets = []

    for islitlet in list_slitlets:

        # define Slitlet2dArc object
        slt = Slitlet2dArc(
            islitlet=islitlet,
            params=params,
            parmodel=parmodel,
            csu_conf=csu_conf,
            ymargin_bb=args_ymargin_bb,
            debugplot=debugplot
        )

        # extract 2D image corresponding to the selected slitlet, clipping
        # the image beyond the unrectified slitlet (in order to isolate
        # the arc lines of the current slitlet; otherwise there are problems
        # with arc lines from neighbour slitlets)
        image2d_tmp = select_unrectified_slitlet(
            image2d=image2d,
            islitlet=islitlet,
            csu_bar_slit_center=csu_conf.csu_bar_slit_center(islitlet),
            params=params,
            parmodel=parmodel,
            maskonly=False
        )
        slitlet2d = slt.extract_slitlet2d(image2d_tmp)

        # subtract smooth background computed as follows:
        # - median collapsed spectrum of the whole slitlet2d
        # - independent median filtering of the previous spectrum in the
        #   two halves in the spectral direction
        if args_remove_sp_background:
            spmedian = np.median(slitlet2d, axis=0)
            naxis1_tmp = spmedian.shape[0]
            jmidpoint = naxis1_tmp // 2
            sp1 = medfilt(spmedian[:jmidpoint], [201])
            sp2 = medfilt(spmedian[jmidpoint:], [201])
            spbackground = np.concatenate((sp1, sp2))
            slitlet2d -= spbackground

        # locate unknown arc lines
        slt.locate_unknown_arc_lines(
            slitlet2d=slitlet2d,
            times_sigma_threshold=args_times_sigma_threshold)

        # continue working with current slitlet only if arc lines have
        # been detected
        if slt.list_arc_lines is not None:

            # compute intersections between spectrum trails and arc lines
            slt.xy_spectrail_arc_intersections(slitlet2d=slitlet2d)

            # compute rectification transformation
            slt.estimate_tt_to_rectify(order=args_order_fmap,
                                       slitlet2d=slitlet2d)

            # rectify image
            slitlet2d_rect = slt.rectify(slitlet2d,
                                         resampling=2,
                                         transformation=1)

            # median spectrum and line peaks from rectified image
            sp_median, fxpeaks = slt.median_spectrum_from_rectified_image(
                slitlet2d_rect,
                sigma_gaussian_filtering=args_sigma_gaussian_filtering,
                nwinwidth_initial=5,
                nwinwidth_refined=5,
                times_sigma_threshold=5,
                npix_avoid_border=6,
                nbrightlines=nbrightlines
            )

            image2d_55sp[islitlet - 1, :] = sp_median

            # determine expected wavelength limits prior to the wavelength
            # calibration
            csu_bar_slit_center = csu_conf.csu_bar_slit_center(islitlet)
            crval1_linear = poly_crval1_linear(csu_bar_slit_center)
            cdelt1_linear = poly_cdelt1_linear(csu_bar_slit_center)
            expected_wvmin = crval1_linear - args_margin_npix * cdelt1_linear
            naxis1_linear = sp_median.shape[0]
            crvaln_linear = crval1_linear + (naxis1_linear - 1) * cdelt1_linear
            expected_wvmax = crvaln_linear + args_margin_npix * cdelt1_linear

            # clip initial master arc line list with bright lines to expected
            # wavelength range
            lok1 = expected_wvmin <= wv_master
            lok2 = wv_master <= expected_wvmax
            lok = lok1 * lok2
            wv_master_eff = wv_master[lok]

            # perform initial wavelength calibration
            solution_wv = wvcal_spectrum(
                sp=sp_median,
                fxpeaks=fxpeaks,
                poly_degree_wfit=args_poldeg_initial,
                wv_master=wv_master_eff,
                wv_ini_search=expected_wvmin,
                wv_end_search=expected_wvmax,
                geometry=args_geometry,
                debugplot=slt.debugplot
            )
            # store initial wavelength calibration polynomial in current
            # slitlet instance
            slt.wpoly = np.polynomial.Polynomial(solution_wv.coeff)
            pause_debugplot(debugplot)

            # clip initial master arc line list with all the lines to expected
            # wavelength range
            lok1 = expected_wvmin <= wv_master_all
            lok2 = wv_master_all <= expected_wvmax
            lok = lok1 * lok2
            wv_master_all_eff = wv_master_all[lok]

            # refine wavelength calibration
            if args_poldeg_refined > 0:
                plottitle = '[slitlet#{}, refined]'.format(islitlet)
                poly_refined, yres_summary = refine_arccalibration(
                    sp=sp_median,
                    poly_initial=slt.wpoly,
                    wv_master=wv_master_all_eff,
                    poldeg=args_poldeg_refined,
                    ntimes_match_wv=1,
                    interactive=args_interactive,
                    threshold=args_threshold_wv,
                    plottitle=plottitle,
                    ylogscale=args_ylogscale,
                    geometry=args_geometry,
                    pdf=args_pdf,
                    debugplot=slt.debugplot
                )
                # store refined wavelength calibration polynomial in current
                # slitlet instance
                slt.wpoly = poly_refined

            # compute approximate linear values for CRVAL1 and CDELT1
            naxis1_linear = sp_median.shape[0]
            crmin1_linear = slt.wpoly(1)
            crmax1_linear = slt.wpoly(naxis1_linear)
            slt.crval1_linear = crmin1_linear
            slt.cdelt1_linear = \
                (crmax1_linear - crmin1_linear) / (naxis1_linear - 1)

            # check that the trimming of wv_master and wv_master_all has
            # preserved the wavelength range [crmin1_linear, crmax1_linear]
            if crmin1_linear < expected_wvmin:
                print(">>> islitlet: ", islitlet)
                print(">>> expected_wvmin:", expected_wvmin)
                print(">>> crmin1_linear.:", crmin1_linear)
                print(">>> WARNING: Unexpected crmin1_linear < expected_wvmin")
            if crmax1_linear > expected_wvmax:
                print(">>> islitlet: ", islitlet)
                print(">>> expected_wvmax:", expected_wvmax)
                print(">>> crmax1_linear.:", crmax1_linear)
                print(">>> WARNING: Unexpected crmax1_linear > expected_wvmin")

            cout = '.'

        else:

            cout = 'x'

        # store current slitlet in list of measured slitlets
        measured_slitlets.append(slt)

        if debugplot == 0:
            if islitlet % 10 == 0:
                if cout != 'x':
                    cout = str(islitlet // 10)
            sys.stdout.write(cout)
            # print(slt)
            if islitlet == list_slitlets[-1]:
                sys.stdout.write('\n')
            sys.stdout.flush()
        else:
            pause_debugplot(debugplot)

    # ---

    # ToDo: return out_55sp ???
    '''
    # Save image with collapsed spectra employed to determine the
    # wavelength calibration
    if args.out_55sp is not None:
        save_ndarray_to_fits(
            array=image2d_55sp,
            file_name=args.out_55sp,
            cast_to_float=True,
            overwrite=True
        )
    '''

    # ---

    # Generate structure to store intermediate results
    outdict = {}
    outdict['instrument'] = 'EMIR'
    outdict['meta-info'] = {}
    outdict['meta-info']['creation_date'] = datetime.now().isoformat()
    outdict['meta-info']['description'] = \
        'computation of rectification and wavelength calibration polynomial ' \
        'coefficients for a particular CSU configuration'
    outdict['meta-info']['recipe_name'] = 'undefined'
    outdict['meta-info']['origin'] = {}
    outdict['meta-info']['origin']['bound_param_uuid'] = \
        bound_param.uuid
    outdict['meta-info']['origin']['arc_image_uuid'] = 'undefined'
    outdict['tags'] = {}
    outdict['tags']['grism'] = grism_name
    outdict['tags']['filter'] = filter_name
    outdict['tags']['islitlet_min'] = islitlet_min
    outdict['tags']['islitlet_max'] = islitlet_max
    outdict['dtu_configuration'] = dtu_conf.outdict()
    outdict['uuid'] = str(uuid4())
    outdict['contents'] = {}

    for slt in measured_slitlets:

        islitlet = slt.islitlet

        # avoid error when creating a python list of coefficients from
        # numpy polynomials when the polynomials do not exist (note that the
        # JSON format doesn't handle numpy arrays and such arrays must be
        # transformed into native python lists)
        if slt.wpoly is None:
            wpoly_coeff = None
        else:
            wpoly_coeff = slt.wpoly.coef.tolist()
        if slt.wpoly_longslit_model is None:
            wpoly_coeff_longslit_model = None
        else:
            wpoly_coeff_longslit_model = slt.wpoly_longslit_model.coef.tolist()

        # avoid similar error when creating a python list of coefficients
        # when the numpy array does not exist; note that this problem
        # does not happen with tt?_aij_longslit_model and
        # tt?_bij_longslit_model because the latter have already been created
        # as native python lists
        if slt.ttd_aij is None:
            ttd_aij = None
        else:
            ttd_aij = slt.ttd_aij.tolist()
        if slt.ttd_bij is None:
            ttd_bij = None
        else:
            ttd_bij = slt.ttd_bij.tolist()
        if slt.tti_aij is None:
            tti_aij = None
        else:
            tti_aij = slt.tti_aij.tolist()
        if slt.tti_bij is None:
            tti_bij = None
        else:
            tti_bij = slt.tti_bij.tolist()

        # creating temporary dictionary with the information corresponding to
        # the current slitlett that will be saved in the JSON file
        tmp_dict = {
            'csu_bar_left': slt.csu_bar_left,
            'csu_bar_right': slt.csu_bar_right,
            'csu_bar_slit_center': slt.csu_bar_slit_center,
            'csu_bar_slit_width': slt.csu_bar_slit_width,
            'x0_reference': slt.x0_reference,
            'y0_reference_lower': slt.y0_reference_lower,
            'y0_reference_middle': slt.y0_reference_middle,
            'y0_reference_upper': slt.y0_reference_upper,
            'y0_frontier_lower': slt.y0_frontier_lower,
            'y0_frontier_upper': slt.y0_frontier_upper,
            'ymargin_bb': slt.ymargin_bb,
            'bb_nc1_orig': slt.bb_nc1_orig,
            'bb_nc2_orig': slt.bb_nc2_orig,
            'bb_ns1_orig': slt.bb_ns1_orig,
            'bb_ns2_orig': slt.bb_ns2_orig,
            'spectrail': {
                'poly_coef_lower':
                    slt.list_spectrails[
                        slt.i_lower_spectrail].poly_funct.coef.tolist(),
                'poly_coef_middle':
                    slt.list_spectrails[
                        slt.i_middle_spectrail].poly_funct.coef.tolist(),
                'poly_coef_upper':
                    slt.list_spectrails[
                        slt.i_upper_spectrail].poly_funct.coef.tolist(),
            },
            'frontier': {
                'poly_coef_lower':
                    slt.list_frontiers[0].poly_funct.coef.tolist(),
                'poly_coef_upper':
                    slt.list_frontiers[1].poly_funct.coef.tolist(),
            },
            'ttd_order': slt.ttd_order,
            'ttd_aij': ttd_aij,
            'ttd_bij': ttd_bij,
            'tti_aij': tti_aij,
            'tti_bij': tti_bij,
            'ttd_order_longslit_model': slt.ttd_order_longslit_model,
            'ttd_aij_longslit_model': slt.ttd_aij_longslit_model,
            'ttd_bij_longslit_model': slt.ttd_bij_longslit_model,
            'tti_aij_longslit_model': slt.tti_aij_longslit_model,
            'tti_bij_longslit_model': slt.tti_bij_longslit_model,
            'wpoly_coeff': wpoly_coeff,
            'wpoly_coeff_longslit_model': wpoly_coeff_longslit_model,
            'crval1_linear': slt.crval1_linear,
            'cdelt1_linear': slt.cdelt1_linear
        }

        slitlet_label = "slitlet" + str(islitlet).zfill(2)
        outdict['contents'][slitlet_label] = tmp_dict

    # ---

    # OBSOLETE
    '''
    # save JSON file needed to compute the MOS model
    with open(args.out_json.name, 'w') as fstream:
        json.dump(outdict, fstream, indent=2, sort_keys=True)
        print('>>> Saving file ' + args.out_json.name)
    '''

    # ---

    # Create object of type RectWaveCoeff with coefficients for
    # rectification and wavelength calibration
    rectwv_coeff = RectWaveCoeff(instrument='EMIR')
    rectwv_coeff.quality_control = numina.types.qc.QC.GOOD
    rectwv_coeff.tags['grism'] = grism_name
    rectwv_coeff.tags['filter'] = filter_name
    rectwv_coeff.meta_info['origin']['bound_param'] = \
        'uuid' + bound_param.uuid
    rectwv_coeff.meta_info['dtu_configuration'] = outdict['dtu_configuration']
    rectwv_coeff.total_slitlets = EMIR_NBARS
    for i in range(EMIR_NBARS):
        islitlet = i + 1
        dumdict = {'islitlet': islitlet}
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        if cslitlet in outdict['contents']:
            dumdict.update(outdict['contents'][cslitlet])
        else:
            dumdict.update({
                'ttd_order': 0
            })
            rectwv_coeff.missing_slitlets.append(islitlet)
        rectwv_coeff.contents.append(dumdict)
    # debugging __getstate__ and __setstate__
    # rectwv_coeff.writeto(args.out_json.name)
    # print('>>> Saving file ' + args.out_json.name)
    # check_setstate_getstate(rectwv_coeff, args.out_json.name)
    logger.info('Generating RectWaveCoeff object with uuid=' +
                rectwv_coeff.uuid)

    return rectwv_coeff


def rectwv_coeff_from_mos_library(reduced_image, master_rectwv,
                                           ignore_DTUconf=True):
    """Evaluate rect.+wavecal. coefficients from MOS library

    Parameters
    ----------
    reduced_image : HDUList object
        Image with preliminary basic reduction: bpm, bias, dark and
        flatfield.
    master_rectwv : MasterRectWave instance
        Rectification and Wavelength Calibrartion Library product.
        Contains the library of polynomial coefficients necessary
        to generate an instance of RectWaveCoeff with the rectification
        and wavelength calibration coefficients for the particular
        CSU configuration.
    ignore_DTUconf : bool
        If True, ignore differences in DTU configuration.

    Returns
    -------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.

    """

    logger = logging.getLogger(__name__)

    # header
    header = reduced_image[0].header

    # read the CSU configuration from the image header
    csu_conf = CsuConfiguration.define_from_header(header)
    logger.debug(csu_conf)

    # read the DTU configuration from the image header
    dtu_conf = DtuConfiguration.define_from_header(header)
    logger.debug(dtu_conf)

    # retrieve DTU configuration from MasterRectWave object
    dtu_conf_calib = DtuConfiguration.define_from_dictionary(
        master_rectwv.meta_info['dtu_configuration']
    )
    # check that the DTU configuration employed to obtain the calibration
    # corresponds to the DTU configuration in the input FITS file
    if dtu_conf != dtu_conf_calib:
        logger.info('DTU configuration from image header:')
        logger.info(dtu_conf)
        logger.info('DTU configuration from master calibration:')
        logger.info(dtu_conf_calib)
        if ignore_DTUconf:
            logger.warning('DTU configuration differences found!')
        else:
            raise ValueError("DTU configurations do not match!")
    else:
        logger.info('DTU configuration match!')

    # check grism and filter
    filter_name = header['filter']
    logger.debug('Filter: ' + filter_name)
    if filter_name != master_rectwv.tags['filter']:
        raise ValueError('Filter name does not match!')
    grism_name = header['grism']
    logger.debug('Grism: ' + grism_name)
    if grism_name != master_rectwv.tags['grism']:
        raise ValueError('Grism name does not match!')

    # valid slitlet numbers
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in master_rectwv.missing_slitlets:
        list_valid_islitlets.remove(idel)
    logger.debug('valid slitlet numbers: ' + str(list_valid_islitlets))

    # initialize intermediate dictionary with relevant information
    # (note: this dictionary corresponds to an old structure employed to
    # store the information in a JSON file; this is no longer necessary,
    # but here we reuse that dictionary for convenience)
    outdict = {}
    outdict['instrument'] = 'EMIR'
    outdict['meta-info'] = {}
    outdict['meta-info']['creation_date'] = datetime.now().isoformat()
    outdict['meta-info']['description'] = \
        'computation of rectification and wavelength calibration polynomial ' \
        'coefficients for a particular CSU configuration from a MOS model '
    outdict['meta-info']['recipe_name'] = 'undefined'
    outdict['meta-info']['origin'] = {}
    outdict['meta-info']['origin']['fits_frame_uuid'] = 'TBD'
    outdict['meta-info']['origin']['rect_wpoly_mos_uuid'] = \
        master_rectwv.uuid
    outdict['meta-info']['origin']['fitted_boundary_param_uuid'] = \
        master_rectwv.meta_info['origin']['bound_param']
    outdict['tags'] = {}
    outdict['tags']['grism'] = grism_name
    outdict['tags']['filter'] = filter_name
    outdict['dtu_configuration'] = dtu_conf.outdict()
    outdict['uuid'] = str(uuid4())
    outdict['contents'] = {}

    # compute rectification and wavelength calibration coefficients for each
    # slitlet according to its csu_bar_slit_center value
    for islitlet in list_valid_islitlets:
        cslitlet = 'slitlet' + str(islitlet).zfill(2)

        # csu_bar_slit_center of current slitlet in initial FITS image
        csu_bar_slit_center = csu_conf.csu_bar_slit_center(islitlet)

        # input data structure
        tmpdict = master_rectwv.contents[islitlet - 1]
        list_csu_bar_slit_center = tmpdict['list_csu_bar_slit_center']

        # check extrapolations
        if csu_bar_slit_center < min(list_csu_bar_slit_center):
            logger.warning('extrapolating table with ' + cslitlet)
            logger.warning('minimum tabulated value: ' +
                           str(min(list_csu_bar_slit_center)))
            logger.warning('sought value...........: ' +
                           str(csu_bar_slit_center))
        if csu_bar_slit_center > max(list_csu_bar_slit_center):
            logger.warning('extrapolating table with ' + cslitlet)
            logger.warning('maximum tabulated value: ' +
                           str(max(list_csu_bar_slit_center)))
            logger.warning('sought value...........: ' +
                           str(csu_bar_slit_center))

        # rectification coefficients
        ttd_order = tmpdict['ttd_order']
        ncoef = ncoef_fmap(ttd_order)
        outdict['contents'][cslitlet] = {}
        outdict['contents'][cslitlet]['ttd_order'] = ttd_order
        outdict['contents'][cslitlet]['ttd_order_longslit_model'] = None
        for keycoef in ['ttd_aij', 'ttd_bij', 'tti_aij', 'tti_bij']:
            coef_out = []
            for icoef in range(ncoef):
                ccoef = str(icoef).zfill(2)
                list_cij = tmpdict['list_' + keycoef + '_' + ccoef]
                funinterp_coef = interp1d(list_csu_bar_slit_center,
                                          list_cij,
                                          kind='linear',
                                          fill_value='extrapolate')
                # note: funinterp_coef expects a numpy array
                dum = funinterp_coef([csu_bar_slit_center])
                coef_out.append(dum[0])
            outdict['contents'][cslitlet][keycoef] = coef_out
            outdict['contents'][cslitlet][keycoef + '_longslit_model'] = None

        # wavelength calibration coefficients
        ncoef = tmpdict['wpoly_degree'] + 1
        wpoly_coeff = []
        for icoef in range(ncoef):
            ccoef = str(icoef).zfill(2)
            list_cij = tmpdict['list_wpoly_coeff_' + ccoef]
            funinterp_coef = interp1d(list_csu_bar_slit_center,
                                      list_cij,
                                      kind='linear',
                                      fill_value='extrapolate')
            # note: funinterp_coef expects a numpy array
            dum = funinterp_coef([csu_bar_slit_center])
            wpoly_coeff.append(dum[0])
        outdict['contents'][cslitlet]['wpoly_coeff'] = wpoly_coeff
        outdict['contents'][cslitlet]['wpoly_coeff_longslit_model'] = None

        # update cdelt1_linear and crval1_linear
        wpoly_function = np.polynomial.Polynomial(wpoly_coeff)
        crmin1_linear = wpoly_function(1)
        crmax1_linear = wpoly_function(EMIR_NAXIS1)
        cdelt1_linear = (crmax1_linear - crmin1_linear) / (EMIR_NAXIS1 - 1)
        crval1_linear = crmin1_linear
        outdict['contents'][cslitlet]['crval1_linear'] = crval1_linear
        outdict['contents'][cslitlet]['cdelt1_linear'] = cdelt1_linear

        # update CSU keywords
        outdict['contents'][cslitlet]['csu_bar_left'] = \
            csu_conf.csu_bar_left(islitlet)
        outdict['contents'][cslitlet]['csu_bar_right'] = \
            csu_conf.csu_bar_right(islitlet)
        outdict['contents'][cslitlet]['csu_bar_slit_center'] = \
            csu_conf.csu_bar_slit_center(islitlet)
        outdict['contents'][cslitlet]['csu_bar_slit_width'] = \
            csu_conf.csu_bar_slit_width(islitlet)

    # for each slitlet compute spectrum trails and frontiers using the
    # fitted boundary parameters
    fitted_bound_param_json = {
        'contents': master_rectwv.meta_info['refined_boundary_model']
    }
    parmodel = fitted_bound_param_json['contents']['parmodel']
    fitted_bound_param_json.update({'meta_info': {'parmodel': parmodel}})
    params = bound_params_from_dict(fitted_bound_param_json)
    logger.debug('Fitted boundary parameters:')
    logger.debug(params.pretty_print())
    for islitlet in list_valid_islitlets:
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        # csu_bar_slit_center of current slitlet in initial FITS image
        csu_bar_slit_center = csu_conf.csu_bar_slit_center(islitlet)
        # compute and store x0_reference value
        x0_reference = float(EMIR_NAXIS1) / 2.0 + 0.5
        outdict['contents'][cslitlet]['x0_reference'] = x0_reference
        # compute spectrum trails (lower, middle and upper)
        list_spectrails = expected_distorted_boundaries(
            islitlet, csu_bar_slit_center,
            [0, 0.5, 1], params, parmodel,
            numpts=101, deg=5, debugplot=0
        )
        # store spectrails in output JSON file
        outdict['contents'][cslitlet]['spectrail'] = {}
        for idum, cdum in zip(range(3), ['lower', 'middle', 'upper']):
            outdict['contents'][cslitlet]['spectrail']['poly_coef_' + cdum] = \
                list_spectrails[idum].poly_funct.coef.tolist()
            outdict['contents'][cslitlet]['y0_reference_' + cdum] = \
                list_spectrails[idum].poly_funct(x0_reference)
        # compute frontiers (lower, upper)
        list_frontiers = expected_distorted_frontiers(
            islitlet, csu_bar_slit_center,
            params, parmodel,
            numpts=101, deg=5, debugplot=0
        )
        # store frontiers in output JSON
        outdict['contents'][cslitlet]['frontier'] = {}
        for idum, cdum in zip(range(2), ['lower', 'upper']):
            outdict['contents'][cslitlet]['frontier']['poly_coef_' + cdum] = \
                list_frontiers[idum].poly_funct.coef.tolist()
            outdict['contents'][cslitlet]['y0_frontier_' + cdum] = \
                list_frontiers[idum].poly_funct(x0_reference)

    # store bounding box parameters for each slitlet
    xdum = np.linspace(1, EMIR_NAXIS1, num=EMIR_NAXIS1)
    for islitlet in list_valid_islitlets:
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        # parameters already available in the input JSON file
        for par in ['bb_nc1_orig', 'bb_nc2_orig', 'ymargin_bb']:
            outdict['contents'][cslitlet][par] = \
                master_rectwv.contents[islitlet - 1][par]
        # estimate bb_ns1_orig and bb_ns2_orig using the already computed
        # frontiers and the value of ymargin_bb, following the same approach
        # employed in the script rect_wpoly_from_longslit; see
        # Slitlet2dLongSlitArc.__init__()
        poly_lower_frontier = np.polynomial.Polynomial(
            outdict['contents'][cslitlet]['frontier']['poly_coef_lower']
        )
        poly_upper_frontier = np.polynomial.Polynomial(
            outdict['contents'][cslitlet]['frontier']['poly_coef_upper']
        )
        ylower = poly_lower_frontier(xdum)
        yupper = poly_upper_frontier(xdum)
        ymargin_bb = master_rectwv.contents[islitlet - 1]['ymargin_bb']
        bb_ns1_orig = int(ylower.min() + 0.5) - ymargin_bb
        if bb_ns1_orig < 1:
            bb_ns1_orig = 1
        bb_ns2_orig = int(yupper.max() + 0.5) + ymargin_bb
        if bb_ns2_orig > EMIR_NAXIS2:
            bb_ns2_orig = EMIR_NAXIS2
        outdict['contents'][cslitlet]['bb_ns1_orig'] = bb_ns1_orig
        outdict['contents'][cslitlet]['bb_ns2_orig'] = bb_ns2_orig

    # ---

    # Create object of type RectWaveCoeff with coefficients for
    # rectification and wavelength calibration
    rectwv_coeff = RectWaveCoeff(instrument='EMIR')
    rectwv_coeff.quality_control = numina.types.qc.QC.GOOD
    rectwv_coeff.tags['grism'] = grism_name
    rectwv_coeff.tags['filter'] = filter_name
    rectwv_coeff.meta_info['origin']['bound_param'] = \
        master_rectwv.meta_info['origin']['bound_param']
    rectwv_coeff.meta_info['origin']['master_rectwv'] = \
        'uuid' + master_rectwv.uuid
    rectwv_coeff.meta_info['dtu_configuration'] = outdict['dtu_configuration']
    rectwv_coeff.total_slitlets = EMIR_NBARS
    for i in range(EMIR_NBARS):
        islitlet = i + 1
        dumdict = {'islitlet': islitlet}
        cslitlet = 'slitlet' + str(islitlet).zfill(2)
        if cslitlet in outdict['contents']:
            dumdict.update(outdict['contents'][cslitlet])
        else:
            dumdict.update({
                'ttd_order': 0
            })
            rectwv_coeff.missing_slitlets.append(islitlet)
        rectwv_coeff.contents.append(dumdict)
    # debugging __getstate__ and __setstate__
    # rectwv_coeff.writeto(args.out_rect_wpoly.name)
    # print('>>> Saving file ' + args.out_rect_wpoly.name)
    # check_setstate_getstate(rectwv_coeff, args.out_rect_wpoly.name)
    logger.info('Generating RectWaveCoeff object with uuid=' +
                rectwv_coeff.uuid)

    return rectwv_coeff


def apply_rectwv_coeff(reduced_image, rectwv_coeff,
                       resampling=2, ignore_DTUconf=True,
                       debugplot=0):
    """Compute rectification and wavelength calibration coefficients.

        Parameters
        ----------
        reduced_image : HDUList object
            Image with preliminary basic reduction: bpm, bias, dark and
            flatfield.
        rectwv_coeff : RectWaveCoeff instance
            Rectification and wavelength calibration coefficients for the
            particular CSU configuration.
        resampling : int
            1: nearest neighbour, 2: flux preserving interpolation.
        ignore_DTUconf : bool
            If True, ignore differences in DTU configuration.
        debugplot : int
            Debugging level for messages and plots. For details see
            'numina.array.display.pause_debugplot.py'.

        Returns
        -------
        rectwv_image : HDUList object
            Rectified and wavelength calibrated image.

        """

    logger = logging.getLogger(__name__)

    # header and data array
    header = reduced_image[0].header
    image2d = reduced_image[0].data

    # check grism and filter
    filter_name = header['filter']
    logger.debug('Filter: ' + filter_name)
    if filter_name != rectwv_coeff.tags['filter']:
        raise ValueError('Filter name does not match!')
    grism_name = header['grism']
    logger.debug('Grism: ' + grism_name)
    if grism_name != rectwv_coeff.tags['grism']:
        raise ValueError('Grism name does not match!')

    # read the DTU configuration from the image header
    dtu_conf = DtuConfiguration.define_from_header(header)
    logger.debug(dtu_conf)

    # retrieve DTU configuration from RectWaveCoeff object
    dtu_conf_calib = DtuConfiguration.define_from_dictionary(
        rectwv_coeff.meta_info['dtu_configuration']
    )
    # check that the DTU configuration employed to obtain the calibration
    # corresponds to the DTU configuration in the input FITS file
    if dtu_conf != dtu_conf_calib:
        logger.info('DTU configuration from image header:')
        logger.info(dtu_conf)
        logger.info('DTU configuration from master calibration:')
        logger.info(dtu_conf_calib)
        if ignore_DTUconf:
            logger.warning('DTU configuration differences found!')
        else:
            raise ValueError("DTU configurations do not match!")
    else:
        logger.info('DTU configuration match!')

    # valid slitlet numbers
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in rectwv_coeff.missing_slitlets:
        list_valid_islitlets.remove(idel)
    logger.debug('Valid slitlet numbers: ' + str(list_valid_islitlets))

    # ---

    # relevant wavelength calibration parameters for rectified and wavelength
    # calibrated image
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    crpix1_enlarged = wv_parameters['crpix1_enlarged']
    crval1_enlarged = wv_parameters['crval1_enlarged']
    cdelt1_enlarged = wv_parameters['cdelt1_enlarged']
    naxis1_enlarged = wv_parameters['naxis1_enlarged']

    # initialize rectified and wavelength calibrated image
    image2d_rectwv = np.zeros((EMIR_NAXIS2, naxis1_enlarged))

    # main loop
    logger.info('Computing rectification and wavelength calibration')
    for islitlet in list_valid_islitlets:
        if abs(debugplot) >= 10:
            if islitlet == list_valid_islitlets[0]:
                print(time.ctime())
            islitlet_progress(islitlet, EMIR_NBARS)
            if islitlet == list_valid_islitlets[-1]:
                print(' ')
                print(time.ctime())

        # define Slitlet2D object
        slt = Slitlet2D(islitlet=islitlet,
                        rectwv_coeff=rectwv_coeff,
                        debugplot=debugplot)

        # extract (distorted) slitlet from the initial image
        slitlet2d = slt.extract_slitlet2d(image2d)

        # rectify slitlet
        slitlet2d_rect = slt.rectify(slitlet2d, resampling=resampling)

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
        image2d_rectwv[i1:i2, :] = slitlet2d_rect_wv[ii1:ii2, :]

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

    logger.info('Updating image header')
    # update wavelength calibration in FITS header
    header.remove('crval1')
    header.remove('crpix1')
    header.remove('crval2')
    header.remove('crpix2')
    header['crpix1'] = (crpix1_enlarged, 'reference pixel')
    header['crval1'] = (crval1_enlarged, 'central wavelength at crpix1')
    header['cdelt1'] = (cdelt1_enlarged, 'linear dispersion (Angstrom/pixel)')
    header['cunit1'] = ('Angstrom', 'units along axis1')
    header['ctype1'] = 'WAVELENGTH'
    header['crpix2'] = (0.0, 'reference pixel')
    header['crval2'] = (0.0, 'central value at crpix2')
    header['cdelt2'] = (1.0, 'increment')
    header['ctype2'] = 'PIXEL'
    header['cunit2'] = ('Pixel', 'units along axis2')
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
    header['history'] = 'Boundary parameters uuid:' + \
                        rectwv_coeff.meta_info['origin']['bound_param'][4:]
    if 'master_rectwv' in rectwv_coeff.meta_info['origin']:
        header['history'] = \
            'MasterRectWave uuid:' + \
            rectwv_coeff.meta_info['origin']['master_rectwv'][4:]
    header['history'] = 'RectWaveCoeff uuid:' + rectwv_coeff.uuid
    header['history'] = 'Rectification and wavelength calibration time ' \
                        + datetime.now().isoformat()

    logger.info('Generating rectified and wavelength calibrated image')

    rectwv_image = fits.PrimaryHDU(data=image2d_rectwv, header=header)

    return rectwv_image
