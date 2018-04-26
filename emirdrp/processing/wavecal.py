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
import numpy as np
from scipy.interpolate import interp1d
import sys
import time
from uuid import uuid4

from numina.array.display.ximshow import ximshow
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.distortion import ncoef_fmap
from numina.array.distortion import order_fmap
from numina.array.distortion import rectify2d
from numina.array.wavecalib.resample import resample_image2d_flux
import numina.types.qc

from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.instrument.dtu_configuration import DtuConfiguration
from emirdrp.products import RectWaveCoeff
from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.fit_boundaries import expected_distorted_boundaries
from emirdrp.tools.fit_boundaries import expected_distorted_frontiers

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_VALID_FILTERS
from emirdrp.core import EMIR_VALID_GRISMS


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


def nscan_minmax_frontiers(y0_frontier_lower, y0_frontier_upper,
                           resize=False):
    """Compute valid scan range for provided y0_frontier values.

    Parameters
    ----------
    y0_frontier_lower : float
        Ordinate of the lower frontier.
    y0_frontier_upper : float
        Ordinate of the upper frontier.
    resize : bool
        If True, when the limits are beyond the expected values
        [1,EMIR_NAXIS2], the values are truncated.

    Returns
    -------
    nscan_min : int
        Minimum useful scan for the image.
    nscan_max : int
        Maximum useful scan for the image.

    """

    fraction_pixel = y0_frontier_lower - int(y0_frontier_lower)
    if fraction_pixel > 0.0:
        nscan_min = int(y0_frontier_lower) + 1
    else:
        nscan_min = int(y0_frontier_lower)
    if nscan_min < 1:
        if resize:
            nscan_min = 1
        else:
            raise ValueError("nscan_min=" + str(nscan_min) + " is < 1")

    fraction_pixel = y0_frontier_upper - int(y0_frontier_upper)
    if fraction_pixel > 0.0:
        nscan_max = int(y0_frontier_upper)
    else:
        nscan_max = int(y0_frontier_upper) - 1
    if nscan_max > EMIR_NAXIS2:
        if resize:
            nscan_max = EMIR_NAXIS2
        else:
            raise ValueError("nscan_max=" + str(nscan_max) +
                             " is > NAXIS2_EMIR=" + str(EMIR_NAXIS2))

    return nscan_min, nscan_max


def evaluate_rectwv_coeff(reduced_image, master_rectwv,
                          ignore_DTUconf=True):
    """Compute rectification and wavelength calibration coefficients.

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

    logger.info('Generating rectified and wavelength calibrated image')

    rectwv_image = fits.PrimaryHDU(data=image2d_rectwv, header=header)

    return rectwv_image
