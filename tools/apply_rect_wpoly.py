from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
import json
import numpy as np
import sys

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow
from numina.array.distortion import order_fmap
from numina.array.distortion import rectify2d
from numina.array.wavecalib.resample import resample_image2d_flux

from arg_file_is_new import arg_file_is_new
from emirdrp.instrument.dtu_configuration import DtuConfiguration
from nscan_minmax_frontiers import nscan_minmax_frontiers
from rect_wpoly_for_mos import islitlet_progress
from save_ndarray_to_fits import save_ndarray_to_fits
from set_wv_parameters import set_wv_parameters

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


class Slitlet2D(object):
    """Slitlet2D class definition.

    Parameters
    ----------
    islitlet : int
        Slitlet number.
    megadict : dictionary
        Python dictionary storing the JSON input file where the
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

    def __init__(self, islitlet, megadict, debugplot):
        # slitlet number
        self.islitlet = islitlet
        cslitlet = 'slitlet' + str(islitlet).zfill(2)

        tmpcontent = megadict['contents'][cslitlet]

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


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(prog='apply_rect_wpoly')

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file",
                        type=argparse.FileType('r'))
    parser.add_argument("--coef_rect_wpoly", required=True,
                        help="Input JSON file with rectification and "
                             "wavelength calibration coefficients",
                        type=argparse.FileType('r'))
    parser.add_argument("--outfile", required=True,
                        help="Output FITS file with rectified and "
                             "wavelength calibrated image",
                        type=lambda x: arg_file_is_new(parser, x))

    # optional arguments
    parser.add_argument("--ignore_DTUconf",
                        help="Ignore DTU configurations differences between "
                             "transformation and input image",
                        action="store_true")
    parser.add_argument("--outfile_rectified_only",
                        help="Output FITS file with rectified image (not"
                             "wavelength calibrated)",
                        type=lambda x: arg_file_is_new(parser, x))
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting & debugging options"
                             " (default=0)",
                        default=0, type=int,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")
    args = parser.parse_args(args)

    if args.echo:
        print('\033[1m\033[31m% ' + ' '.join(sys.argv) + '\033[0m\n')

    # read calibration structure from JSON file
    rect_wpoly_dict = json.loads(open(args.coef_rect_wpoly.name).read())

    # read FITS image and its corresponding header
    hdulist = fits.open(args.fitsfile)
    header = hdulist[0].header
    image2d = hdulist[0].data
    hdulist.close()

    # protections
    naxis2, naxis1 = image2d.shape
    if naxis1 != header['naxis1'] or naxis2 != header['naxis2']:
        print('>>> NAXIS1:', naxis1)
        print('>>> NAXIS2:', naxis2)
        raise ValueError('Something is wrong with NAXIS1 and/or NAXIS2')
    if abs(args.debugplot) >= 10:
        print('>>> NAXIS1:', naxis1)
        print('>>> NAXIS2:', naxis2)

    # check that the input FITS file grism and filter match
    filter_name = header['filter']
    if filter_name != rect_wpoly_dict['tags']['filter']:
        raise ValueError("Filter name does not match!")
    grism_name = header['grism']
    if grism_name != rect_wpoly_dict['tags']['grism']:
        raise ValueError("Filter name does not match!")
    if abs(args.debugplot) >= 10:
        print('>>> grism.......:', grism_name)
        print('>>> filter......:', filter_name)

    # check that the DTU configurations are compatible
    dtu_conf_fitsfile = DtuConfiguration.define_from_fits(args.fitsfile)
    dtu_conf_jsonfile = DtuConfiguration.define_from_dictionary(
        rect_wpoly_dict['dtu_configuration'])
    if dtu_conf_fitsfile != dtu_conf_jsonfile:
        print('DTU configuration (FITS file):\n\t', dtu_conf_fitsfile)
        print('DTU configuration (JSON file):\n\t', dtu_conf_jsonfile)
        if args.ignore_DTUconf:
            print('WARNING: DTU configuration differences found!')
        else:
            raise ValueError("DTU configurations do not match!")
    else:
        if abs(args.debugplot) >= 10:
            print('>>> DTU Configuration match!')
            print(dtu_conf_fitsfile)

    # read islitlet_min and islitlet_max from input JSON file
    islitlet_min = rect_wpoly_dict['tags']['islitlet_min']
    islitlet_max = rect_wpoly_dict['tags']['islitlet_max']
    if abs(args.debugplot) >= 10:
        print('>>> islitlet_min:', islitlet_min)
        print('>>> islitlet_max:', islitlet_max)

    # ---

    if args.outfile_rectified_only is not None:
        image2d_rectified = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))
        image2d_unrectified = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

        for islitlet in range(islitlet_min, islitlet_max + 1):
            if args.debugplot == 0:
                islitlet_progress(islitlet, islitlet_max)

            # define Slitlet2D object
            slt = Slitlet2D(islitlet=islitlet,
                            megadict=rect_wpoly_dict,
                            debugplot=args.debugplot)

            # minimum and maximum useful scan (pixel in the spatial direction)
            # for the rectified slitlet
            nscan_min, nscan_max = nscan_minmax_frontiers(
                slt.y0_frontier_lower,
                slt.y0_frontier_upper,
                resize=False
            )
            # extract 2D image corresponding to the selected slitlet: note that
            # in this case we are not using select_unrectified_slitlets()
            # because it introduces extra zero pixels in the slitlet frontiers
            slitlet2d = slt.extract_slitlet2d(image2d)

            # rectify image
            slitlet2d_rect = slt.rectify(slitlet2d,
                                         resampling=1)
            slitlet2d_unrect = slt.rectify(slitlet2d_rect,
                                           resampling=1,
                                           inverse=True)

            ii1 = nscan_min - slt.bb_ns1_orig
            ii2 = nscan_max - slt.bb_ns1_orig + 1

            j1 = slt.bb_nc1_orig - 1
            j2 = slt.bb_nc2_orig
            i1 = slt.bb_ns1_orig - 1 + ii1
            i2 = i1 + ii2 - ii1

            image2d_rectified[i1:i2, j1:j2] = slitlet2d_rect[ii1:ii2, :]
            image2d_unrectified[i1:i2, j1:j2] = slitlet2d_unrect[ii1:ii2, :]

        save_ndarray_to_fits(
            array=[image2d_rectified, image2d_unrectified],
            file_name=args.outfile_rectified_only,
            cast_to_float=[True] * 2,
            overwrite=True
        )

    # ---

    # relevant wavelength calibration parameters for rectified and wavelength
    # calibrated image
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    crpix1_enlarged = wv_parameters['crpix1_enlarged']
    crval1_enlarged = wv_parameters['crval1_enlarged']
    cdelt1_enlarged = wv_parameters['cdelt1_enlarged']
    naxis1_enlarged = wv_parameters['naxis1_enlarged']

    # initialize rectified and wavelength calibrated image
    image2d_rectified_wv = np.zeros((EMIR_NAXIS2, naxis1_enlarged))

    # main loop
    for islitlet in range(islitlet_min, islitlet_max + 1):
        if args.debugplot == 0:
            islitlet_progress(islitlet, islitlet_max)

        # define Slitlet2D object
        slt = Slitlet2D(islitlet=islitlet,
                        megadict=rect_wpoly_dict,
                        debugplot=args.debugplot)

        if abs(args.debugplot) >= 10:
            print(slt)

        # extract (distorted) slitlet from the initial image
        slitlet2d = slt.extract_slitlet2d(image2d)

        # rectify slitlet
        slitlet2d_rect = slt.rectify(slitlet2d, resampling=1)

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
        image2d_rectified_wv[i1:i2, :] = slitlet2d_rect_wv[ii1:ii2, :]

        # include scan range in FITS header
        header['sltmin' + str(islitlet).zfill(2)] = i1
        header['sltmax' + str(islitlet).zfill(2)] = i2 - 1

    # modify upper limit of previous slitlet in case of overlapping:
    # note that the overlapped scans have been overwritten with the
    # information from the current slitlet!
    for islitlet in range(islitlet_min, islitlet_max + 1):
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

    save_ndarray_to_fits(
        array=image2d_rectified_wv,
        file_name=args.outfile,
        main_header=header,
        overwrite=True
    )
    print('>>> Saving file ' + args.outfile.name)


if __name__ == "__main__":
    main()
