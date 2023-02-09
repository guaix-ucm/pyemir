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

"""Definition of Slitlet2D class"""

from __future__ import division, print_function

import numpy as np

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow
from numina.array.distortion import order_fmap
from numina.array.distortion import rectify2d

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NPIXPERSLIT_RECTIFIED


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
    iminslt : int
        Minimum useful scan in output full 2d rectified image.
    imaxslt : int
        Maximum useful scan in output full 2d rectified image.
    jminslt : int
        Minimum useful channel in output full 2d rectified image.
    jmaxslt : int
        Maximum useful channel in output full 2d rectified image.
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

        # define expected frontier ordinates at x0_reference and
        # associated linear transformation
        self.y0_frontier_lower_expected = \
            tmpcontent['y0_frontier_lower_expected']
        self.y0_frontier_upper_expected = \
            tmpcontent['y0_frontier_upper_expected']
        self.corr_yrect_a = tmpcontent['corr_yrect_a']
        self.corr_yrect_b = tmpcontent['corr_yrect_b']
        self.min_row_rectified = tmpcontent['min_row_rectified']
        self.max_row_rectified = tmpcontent['max_row_rectified']

        # define useful scan region in output rectified image
        self.iminslt = (islitlet - 1) * EMIR_NPIXPERSLIT_RECTIFIED + 1
        self.imaxslt = islitlet * EMIR_NPIXPERSLIT_RECTIFIED

        # define useful channel region in output rectified image
        self.jminslt = 0
        self.jmaxslt = 0

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
                 "- islitlet..................: " + \
                 str(self.islitlet) + "\n" + \
                 "- csu_bar_left..............: " + \
                 str(self.csu_bar_left) + "\n" + \
                 "- csu_bar_right.............: " + \
                 str(self.csu_bar_right) + "\n" + \
                 "- csu_bar_slit_center.......: " + \
                 str(self.csu_bar_slit_center) + "\n" + \
                 "- csu_bar_slit_width........: " + \
                 str(self.csu_bar_slit_width) + "\n" + \
                 "- x0_reference..............: " + \
                 str(self.x0_reference) + "\n" + \
                 "- y0_reference_lower........: " + \
                 str(self.y0_reference_lower) + "\n" + \
                 "- y0_reference_middle.......: " + \
                 str(self.y0_reference_middle) + "\n" + \
                 "- y0_reference_upper........: " + \
                 str(self.y0_reference_upper) + "\n" + \
                 "- y0_frontier_lower.........: " + \
                 str(self.y0_frontier_lower) + "\n" + \
                 "- y0_frontier_upper.........: " + \
                 str(self.y0_frontier_upper) + "\n" + \
                 "- y0_frontier_lower_expected: " + \
                 str(self.y0_frontier_lower_expected) + "\n" + \
                 "- y0_frontier_upper_expected: " + \
                 str(self.y0_frontier_upper_expected) + "\n" + \
                 "- corr_yrect_a................: " + \
                 str(self.corr_yrect_a) + "\n" + \
                 "- corr_yrect_b................: " + \
                 str(self.corr_yrect_b) + "\n" + \
                 "- min_row_rectified...........: " + \
                 str(self.min_row_rectified) + "\n" + \
                 "- max_row_rectified...........: " + \
                 str(self.max_row_rectified) + "\n" + \
                 "- iminslt.....................: " + \
                 str(self.iminslt) + "\n" + \
                 "- imaxslt.....................: " + \
                 str(self.imaxslt) + "\n" + \
                 "- jminslt.....................: " + \
                 str(self.jminslt) + "\n" + \
                 "- jmaxslt.....................: " + \
                 str(self.jmaxslt) + "\n" + \
                 "- bb_nc1_orig...............: " + \
                 str(self.bb_nc1_orig) + "\n" + \
                 "- bb_nc2_orig...............: " + \
                 str(self.bb_nc2_orig) + "\n" + \
                 "- bb_ns1_orig...............: " + \
                 str(self.bb_ns1_orig) + "\n" + \
                 "- bb_ns2_orig...............: " + \
                 str(self.bb_ns2_orig) + "\n" + \
                 "- lower spectrail...........:\n\t" + \
                 str(self.list_spectrails[0]) + "\n" + \
                 "- middle spectrail..........:\n\t" + \
                 str(self.list_spectrails[1]) + "\n" + \
                 "- upper spectrail...........:\n\t" + \
                 str(self.list_spectrails[2]) + "\n" + \
                 "- lower frontier............:\n\t" + \
                 str(self.list_frontiers[0]) + "\n" + \
                 "- upper frontier............:\n\t" + \
                 str(self.list_frontiers[1]) + "\n" + \
                 "- ttd_order.................: " + \
                 str(self.ttd_order) + "\n" + \
                 "- ttd_aij...................:\n\t" + \
                 str(self.ttd_aij) + "\n" + \
                 "- ttd_bij...................:\n\t" + \
                 str(self.ttd_bij) + "\n" + \
                 "- tti_aij...................:\n\t" + \
                 str(self.tti_aij) + "\n" + \
                 "- tti_bij...................:\n\t" + \
                 str(self.tti_bij) + "\n" + \
                 "- wpoly.....................:\n\t" + \
                 str(self.wpoly) + "\n" + \
                 "- debugplot.................: " + \
                 str(self.debugplot)

        return output

    def extract_slitlet2d(self, image_2k2k, subtitle=None):
        """Extract slitlet 2d image from image with original EMIR dimensions.

        Parameters
        ----------
        image_2k2k : numpy array
            Original image (dimensions EMIR_NAXIS1 * EMIR_NAXIS2)
        subtitle : string, optional
            Subtitle for plot.

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
            self.ximshow_unrectified(slitlet2d, subtitle)

        # return slitlet image
        return slitlet2d

    def rectify(self, slitlet2d, resampling, inverse=False, subtitle=None):
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
        subtitle : string, optional
            Subtitle for plot.

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
                self.ximshow_unrectified(slitlet2d_rect, subtitle=subtitle)
            else:
                self.ximshow_rectified(slitlet2d_rect, subtitle=subtitle)

        return slitlet2d_rect

    def ximshow_unrectified(self, slitlet2d, subtitle=None):
        """Display unrectified image with spectrails and frontiers.

        Parameters
        ----------
        slitlet2d : numpy array
            Array containing the unrectified slitlet image.
        subtitle : string, optional
            Subtitle for plot.

        """

        title = "Slitlet#" + str(self.islitlet)
        if subtitle is not None:
            title += ' (' + subtitle + ')'
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
        if title is not None:
            ax.set_title(title)
        pause_debugplot(debugplot=self.debugplot, pltshow=True)

    def ximshow_rectified(self, slitlet2d_rect, subtitle=None):
        """Display rectified image with spectrails and frontiers.

        Parameters
        ----------
        slitlet2d_rect : numpy array
            Array containing the rectified slitlet image.
        subtitle : string, optional
            Subtitle for plot.

        """

        title = "Slitlet#" + str(self.islitlet)
        if subtitle is not None:
            title += ' (' + subtitle + ')'
        ax = ximshow(slitlet2d_rect, title=title,
                     first_pixel=(self.bb_nc1_orig, self.bb_ns1_orig),
                     show=False)
        # grid with fitted transformation: spectrum trails
        xx = np.arange(0, self.bb_nc2_orig - self.bb_nc1_orig + 1,
                       dtype=float)
        for spectrail in self.list_spectrails:
            yy0 = self.corr_yrect_a + \
                  self.corr_yrect_b * spectrail(self.x0_reference)
            yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
            ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b")
        for spectrail in self.list_frontiers:
            yy0 = self.corr_yrect_a +\
                  self.corr_yrect_b * spectrail(self.x0_reference)
            yy = np.tile([yy0 - self.bb_ns1_orig], xx.size)
            ax.plot(xx + self.bb_nc1_orig, yy + self.bb_ns1_orig, "b:")
        # show plot
        pause_debugplot(self.debugplot, pltshow=True)
