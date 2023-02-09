#
# Copyright 2018-2023 Universidad Complutense de Madrid
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

"""Definition of Slitlet2dArc class"""

from __future__ import division, print_function

from matplotlib.patches import Rectangle
import numpy as np
from scipy import ndimage
from skimage import restoration

from numina.array.ccd_line import ArcLine
from numina.array.ccd_line import intersection_spectrail_arcline
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.polfit_residuals import \
    polfit_residuals_with_sigma_rejection
from numina.array.display.ximplotxy import ximplotxy
from numina.array.display.ximshow import ximshow
from numina.array.distortion import compute_distortion
from numina.array.distortion import fmap
from numina.array.distortion import rectify2d
from numina.array.wavecalib.peaks_spectrum import find_peaks_spectrum
from numina.array.wavecalib.peaks_spectrum import refine_peaks_spectrum

from emirdrp.tools.fit_boundaries import expected_distorted_boundaries
from emirdrp.tools.fit_boundaries import expected_distorted_frontiers

from .rescale_array_z1z2 import rescale_array_to_z1z2
from .rescale_array_z1z2 import rescale_array_from_z1z2

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NPIXPERSLIT_RECTIFIED


def expected_y0_lower_frontier(islitlet):
    """Expected ordinate of lower frontier in rectified image.

    Parameters
    ----------
    islitlet : int
        Slitlet number.

    Returns
    ----------
    y0_lower : float
        Ordinate value.
    """

    y0_lower =  0.5 + float((islitlet - 1) * EMIR_NPIXPERSLIT_RECTIFIED)
    return y0_lower


def expected_y0_upper_frontier(islitlet):
    """Expected ordinate of upper frontier in rectified image.

    Parameters
    ----------
    islitlet : int
        Slitlet number.

    Returns
    ----------
    y0_upper : float
        Ordinate value.
    """

    y0_upper =  float((islitlet) * EMIR_NPIXPERSLIT_RECTIFIED) + 0.5
    return y0_upper


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
    csu_conf : CsuConfiguration object
        Instance of CsuConfiguration.
    ymargin_bb : int
        Extra number of pixels above and below the enclosing rectangle
        defined by the slitlet frontiers.
    params : :class:`~lmfit.parameter.Parameters` or None
        Parameters to be employed in the prediction of the distorted
        boundaries.
    parmodel : str or None
        Model to be assumed. Allowed values are 'longslit' and
        'multislit'.
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
    y0_reference_lower_expected: float
        Expected Y coordinate corresponding to the lower boundary
        computed at x0_reference in the rectified image.
    y0_reference_middle_expected: float
        Expected Y coordinate corresponding to the middle boundary
        computed at x0_reference in the rectified image.
    y0_reference_upper_expected: float
        Expected Y coordinate corresponding to the upper boundary
        computed at x0_reference in the rectified image.
    y0_frontier_lower_expected: float
        Expected Y coordinate corresponding to the lower frontier
        computed at x0_reference in the rectified image.
    y0_frontier_upper_expected: float
        Expected Y coordinate corresponding to the upper frontier
        computed at x0_reference in the rectified image.
    corr_yrect_a : float
        Intercept of the relation y_expected = a + b * y_measured
        that transforms the measured ordinate into the expected ordinate
        of the rectified image.
    corr_yrect_b : float
        Slope of the relation y_expected = a + b * y_measured
        that transforms the measured ordinate into the expected ordinate
        of the rectified image.
    min_row_rectified : int
        Minimum useful row (starting from zero) of the rectified slitlet.
    max_row_rectified : int
        Maximum useful row (starting from zero) of the rectified slitlet.
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

    def __init__(self, islitlet,
                 csu_conf, ymargin_bb,
                 params=None, parmodel=None,
                 debugplot=0):

        # slitlet number
        self.islitlet = islitlet

        # debugplot
        self.debugplot = debugplot

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

        # define expected frontier ordinates at x0_reference for the rectified
        # image imposing the vertical length of the slitlet to be constant
        # and equal to EMIR_NPIXPERSLIT_RECTIFIED
        self.y0_frontier_lower_expected = expected_y0_lower_frontier(
            islitlet)
        self.y0_frontier_upper_expected = expected_y0_upper_frontier(
            islitlet)

        if params is None:
            return

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

        # compute linear transformation to place the rectified slitlet at
        # the center of the current slitlet bounding box
        xdum1 = self.y0_frontier_lower
        ydum1 = self.y0_frontier_lower_expected
        xdum2 = self.y0_frontier_upper
        ydum2 = self.y0_frontier_upper_expected
        self.corr_yrect_b = (ydum2 - ydum1) / (xdum2 - xdum1)
        self.corr_yrect_a = ydum1 - self.corr_yrect_b * xdum1
        # compute expected location of rectified boundaries
        self.y0_reference_lower_expected = \
            self.corr_yrect_a + self.corr_yrect_b * self.y0_reference_lower
        self.y0_reference_middle_expected = \
            self.corr_yrect_a + self.corr_yrect_b * self.y0_reference_middle
        self.y0_reference_upper_expected = \
            self.corr_yrect_a + self.corr_yrect_b * self.y0_reference_upper
        # shift transformation to center the rectified slitlet within the
        # slitlet bounding box
        ydummid = (ydum1 + ydum2) / 2
        ioffset = int(
            ydummid - (self.bb_ns1_orig + self.bb_ns2_orig) / 2.0)
        self.corr_yrect_a -= ioffset

        # minimum and maximum row in the rectified slitlet encompassing
        # EMIR_NPIXPERSLIT_RECTIFIED pixels
        # a) scan number (in pixels, from 1 to NAXIS2)
        xdum1 = self.corr_yrect_a + \
                self.corr_yrect_b * self.y0_frontier_lower
        xdum2 = self.corr_yrect_a + \
                self.corr_yrect_b * self.y0_frontier_upper
        # b) row number (starting from zero)
        self.min_row_rectified = \
            int((round(xdum1 * 10) + 5) / 10) - self.bb_ns1_orig
        self.max_row_rectified = \
            int((round(xdum2 * 10) - 5) / 10) - self.bb_ns1_orig

        # place holder for still undefined class members
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
                 "- y0_frontier_lower_expected: " + \
                 str(self.y0_frontier_lower_expected) + "\n" + \
                 "- y0_frontier_upper_expected: " + \
                 str(self.y0_frontier_upper_expected) + "\n" + \
                 "- corr_yrect_a................: " + \
                 str(self.corr_yrect_a) + \
                 "- corr_yrect_b................: " + \
                 str(self.corr_yrect_b) + \
                 "- min_row_rectified...........: " + \
                 str(self.min_row_rectified) + \
                 "- max_row_rectified...........: " + \
                 str(self.max_row_rectified) + \
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
        slitlet2d = slitlet2d.astype(float)

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
                # compute expected ordinate y_expected in the rectified
                # image
                y_expected = self.corr_yrect_a + self.corr_yrect_b * \
                             spectrail.y_rectified
                self.y_inter_rect = np.append(self.y_inter_rect, y_expected)
        if abs(self.debugplot) >= 10:
            print('>>> y0_frontier_lower_expected........: ',
                  self.y0_frontier_lower_expected)
            print('>>> y0_frontier_upper_expected........: ',
                  self.y0_frontier_upper_expected)
            print('>>> shifted y0_frontier_upper_expected: ',
                  self.corr_yrect_a +
                  self.corr_yrect_b * self.y0_frontier_lower)
            print('>>> shifted y0_frontier_lower_expected: ',
                  self.corr_yrect_a +
                  self.corr_yrect_b * self.y0_frontier_upper)
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
                           dtype=float)
            for spectrail in self.list_spectrails:
                yy0 = self.corr_yrect_a + \
                      self.corr_yrect_b * spectrail.y_rectified
                yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
                ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b")
                xxx, yyy = fmap(self.ttd_order, self.ttd_aij, self.ttd_bij,
                                xx, yy)
                ax.plot(xxx + self.bb_nc1_orig, yyy + self.bb_ns1_orig, "g")
            # grid with fitted transformation: arc lines
            ylower_line = \
                self.list_spectrails[self.i_lower_spectrail].y_rectified
            ylower_line = self.corr_yrect_a + self.corr_yrect_b * ylower_line
            yupper_line = \
                self.list_spectrails[self.i_upper_spectrail].y_rectified
            yupper_line = self.corr_yrect_a + self.corr_yrect_b * yupper_line
            n_points = int(yupper_line - ylower_line + 0.5) + 1
            yy = np.linspace(ylower_line - self.bb_ns1_orig,
                             yupper_line - self.bb_ns1_orig,
                             num=n_points,
                             dtype=float)
            for arc_line in self.list_arc_lines:
                xline = arc_line.x_rectified - self.bb_nc1_orig
                xx = np.array([xline] * n_points)
                ax.plot(xx + self.bb_nc1_orig,yy + self.bb_ns1_orig, "b" )
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
                           dtype=float)
            for spectrail in self.list_spectrails:
                yy0 = self.corr_yrect_a + \
                      self.corr_yrect_b * spectrail.y_rectified
                yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
                ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b")
            for spectrail in self.list_frontiers:
                yy0 = self.corr_yrect_a + \
                      self.corr_yrect_b * spectrail.y_rectified
                yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
                ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b:")
            # grid with fitted transformation: arc lines
            ylower_line = \
                self.list_spectrails[self.i_lower_spectrail].y_rectified
            ylower_line = self.corr_yrect_a + self.corr_yrect_b * ylower_line
            yupper_line = \
                self.list_spectrails[self.i_upper_spectrail].y_rectified
            yupper_line = self.corr_yrect_a + self.corr_yrect_b * yupper_line
            n_points = int(yupper_line - ylower_line + 0.5) + 1
            yy = np.linspace(ylower_line - self.bb_ns1_orig,
                             yupper_line - self.bb_ns1_orig,
                             num=n_points,
                             dtype=float)
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
