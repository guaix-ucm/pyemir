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
from numpy.polynomial.polynomial import polyval
from scipy import ndimage
import sys

from numina.array.display.ximplotxy import ximplotxy
from numina.array.display.ximshow import ximshow
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.numsplines import AdaptiveLSQUnivariateSpline
from numina.array.wavecalib.apply_integer_offsets import apply_integer_offsets
from numina.tools.arg_file_is_new import arg_file_is_new

from emirdrp.processing.wavecal.slitlet2d import Slitlet2D
from emirdrp.processing.wavecal.set_wv_parameters import set_wv_parameters
from emirdrp.instrument.dtuconf import DtuConf
from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers
from emirdrp.tools.rect_wpoly_for_mos import islitlet_progress
from emirdrp.tools.save_ndarray_to_fits import save_ndarray_to_fits
from emirdrp.products import RectWaveCoeff

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: compute pixel-to-pixel flatfield'
    )

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file (flat ON-OFF)",
                        type=argparse.FileType('rb'))
    parser.add_argument("--rectwv_coeff", required=True,
                        help="Input JSON file with rectification and "
                             "wavelength calibration coefficients",
                        type=argparse.FileType('rt'))
    parser.add_argument("--minimum_slitlet_width_mm", required=True,
                        help="Minimum slitlet width in mm",
                        type=float)
    parser.add_argument("--maximum_slitlet_width_mm", required=True,
                        help="Maximum slitlet width in mm",
                        type=float)
    parser.add_argument("--minimum_fraction", required=True,
                        help="Minimum allowed flatfielding value",
                        type=float, default=0.01)
    parser.add_argument("--minimum_value_in_output",
                        help="Minimum value allowed in output file: pixels "
                             "below this value are set to 1.0 (default=0.01)",
                        type=float, default=0.01)
    parser.add_argument("--maximum_value_in_output",
                        help="Maximum value allowed in output file: pixels "
                             "above this value are set to 1.0 (default=10.0)",
                        type=float, default=10.0)
    parser.add_argument("--nwindow_median",
                        help="Window size to smooth median spectrum in the "
                             "spectral direction",
                        type=int)
    parser.add_argument("--outfile", required=True,
                        help="Output FITS file",
                        type=lambda x: arg_file_is_new(parser, x, mode='wb'))

    # optional arguments
    parser.add_argument("--delta_global_integer_offset_x_pix",
                        help="Delta global integer offset in the X direction "
                             "(default=0)",
                        default=0, type=int)
    parser.add_argument("--delta_global_integer_offset_y_pix",
                        help="Delta global integer offset in the Y direction "
                             "(default=0)",
                        default=0, type=int)
    parser.add_argument("--resampling",
                        help="Resampling method: 1 -> nearest neighbor, "
                             "2 -> linear interpolation (default)",
                        default=2, type=int,
                        choices=(1, 2))
    parser.add_argument("--ignore_DTUconf",
                        help="Ignore DTU configurations differences between "
                             "model and input image",
                        action="store_true")
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

    # This code is obsolete
    raise ValueError('This code is obsolete: use recipe in '
                     'emirdrp/recipes/spec/flatpix2pix.py')

    # read calibration structure from JSON file
    rectwv_coeff = RectWaveCoeff._datatype_load(args.rectwv_coeff.name)

    # modify (when requested) global offsets
    rectwv_coeff.global_integer_offset_x_pix += \
        args.delta_global_integer_offset_x_pix
    rectwv_coeff.global_integer_offset_y_pix += \
        args.delta_global_integer_offset_y_pix

    # read FITS image and its corresponding header
    hdulist = fits.open(args.fitsfile)
    header = hdulist[0].header
    image2d = hdulist[0].data
    hdulist.close()

    # apply global offsets
    image2d = apply_integer_offsets(
        image2d=image2d,
        offx=rectwv_coeff.global_integer_offset_x_pix,
        offy=rectwv_coeff.global_integer_offset_y_pix
    )

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
    if filter_name != rectwv_coeff.tags['filter']:
        raise ValueError("Filter name does not match!")
    grism_name = header['grism']
    if grism_name != rectwv_coeff.tags['grism']:
        raise ValueError("Filter name does not match!")
    if abs(args.debugplot) >= 10:
        print('>>> grism.......:', grism_name)
        print('>>> filter......:', filter_name)

    # check that the DTU configurations are compatible
    with fits.open(args.fitsfile) as hdulist:
        dtu_conf_fitsfile = DtuConf.from_img(hdulist)
    dtu_conf_jsonfile = DtuConf.from_values(
        **rectwv_coeff.meta_info['dtu_configuration']
    )
    if dtu_conf_fitsfile != dtu_conf_jsonfile:
        print('DTU configuration (FITS file):\n\t', dtu_conf_fitsfile)
        print('DTU configuration (JSON file):\n\t', dtu_conf_jsonfile)
        if args.ignore_DTUconf:
            print('WARNING: DTU configuration differences found!')
        else:
            raise ValueError('DTU configurations do not match')
    else:
        if abs(args.debugplot) >= 10:
            print('>>> DTU Configuration match!')
            print(dtu_conf_fitsfile)

    # load CSU configuration
    csu_conf_fitsfile = CsuConfiguration.define_from_fits(args.fitsfile)
    if abs(args.debugplot) >= 10:
        print(csu_conf_fitsfile)

    # valid slitlet numbers
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in rectwv_coeff.missing_slitlets:
        print('-> Removing slitlet (not defined):', idel)
        list_valid_islitlets.remove(idel)
    # filter out slitlets with widths outside valid range
    list_outside_valid_width = []
    for islitlet in list_valid_islitlets:
        slitwidth = csu_conf_fitsfile.csu_bar_slit_width(islitlet)
        if (slitwidth < args.minimum_slitlet_width_mm) or \
                (slitwidth > args.maximum_slitlet_width_mm):
            list_outside_valid_width.append(islitlet)
            print('-> Removing slitlet (invalid width):', islitlet)
    if len(list_outside_valid_width) > 0:
        for idel in list_outside_valid_width:
            list_valid_islitlets.remove(idel)
    print('>>> valid slitlet numbers:\n', list_valid_islitlets)

    # ---

    # compute and store median spectrum (and masked region) for each
    # individual slitlet
    image2d_sp_median = np.zeros((EMIR_NBARS, EMIR_NAXIS1))
    image2d_sp_mask = np.zeros((EMIR_NBARS, EMIR_NAXIS1), dtype=bool)
    for islitlet in list(range(1, EMIR_NBARS + 1)):
        if islitlet in list_valid_islitlets:
            if args.debugplot == 0:
                islitlet_progress(islitlet, EMIR_NBARS, ignore=False)
            # define Slitlet2D object
            slt = Slitlet2D(islitlet=islitlet,
                            rectwv_coeff=rectwv_coeff,
                            debugplot=args.debugplot)

            if abs(args.debugplot) >= 10:
                print(slt)

            # extract (distorted) slitlet from the initial image
            slitlet2d = slt.extract_slitlet2d(
                image_2k2k=image2d,
                subtitle='original image'
            )

            # rectify slitlet
            slitlet2d_rect = slt.rectify(
                slitlet2d=slitlet2d,
                resampling=args.resampling,
                subtitle='original rectified'
            )
            naxis2_slitlet2d, naxis1_slitlet2d = slitlet2d_rect.shape

            if naxis1_slitlet2d != EMIR_NAXIS1:
                print('naxis1_slitlet2d: ', naxis1_slitlet2d)
                print('EMIR_NAXIS1.....: ', EMIR_NAXIS1)
                raise ValueError("Unexpected naxis1_slitlet2d")

            sp_mask = np.zeros(naxis1_slitlet2d, dtype=bool)

            # for grism LR set to zero data beyond useful wavelength range
            if grism_name == 'LR':
                wv_parameters = set_wv_parameters(filter_name, grism_name)
                x_pix = np.arange(1, naxis1_slitlet2d + 1)
                wl_pix = polyval(x_pix, slt.wpoly)
                lremove = wl_pix < wv_parameters['wvmin_useful']
                sp_mask[lremove] = True
                slitlet2d_rect[:, lremove] = 0.0
                lremove = wl_pix > wv_parameters['wvmax_useful']
                slitlet2d_rect[:, lremove] = 0.0
                sp_mask[lremove] = True

            # get useful slitlet region (use boundaries instead of frontiers;
            # note that the nscan_minmax_frontiers() works well independently
            # of using frontiers of boundaries as arguments)
            nscan_min, nscan_max = nscan_minmax_frontiers(
                slt.y0_reference_lower,
                slt.y0_reference_upper,
                resize=False
            )
            ii1 = nscan_min - slt.bb_ns1_orig
            ii2 = nscan_max - slt.bb_ns1_orig + 1

            # median spectrum
            sp_collapsed = np.median(slitlet2d_rect[ii1:(ii2 + 1), :], axis=0)

            # smooth median spectrum along the spectral direction
            sp_median = ndimage.median_filter(
                sp_collapsed,
                args.nwindow_median,
                mode='nearest'
            )

            """
                nremove = 5
                spl = AdaptiveLSQUnivariateSpline(
                    x=xaxis1[nremove:-nremove],
                    y=sp_collapsed[nremove:-nremove],
                    t=11,
                    adaptive=True
                )
                xknots = spl.get_knots()
                yknots = spl(xknots)
                sp_median = spl(xaxis1)

                # compute rms within each knot interval
                nknots = len(xknots)
                rms_array = np.zeros(nknots - 1, dtype=float)
                for iknot in range(nknots - 1):
                    residuals = []
                    for xdum, ydum, yydum in \
                            zip(xaxis1, sp_collapsed, sp_median):
                        if xknots[iknot] <= xdum <= xknots[iknot + 1]:
                            residuals.append(abs(ydum - yydum))
                    if len(residuals) > 5:
                        rms_array[iknot] = np.std(residuals)
                    else:
                        rms_array[iknot] = 0

                # determine in which knot interval falls each pixel
                iknot_array = np.zeros(len(xaxis1), dtype=int)
                for idum, xdum in enumerate(xaxis1):
                    for iknot in range(nknots - 1):
                        if xknots[iknot] <= xdum <= xknots[iknot + 1]:
                            iknot_array[idum] = iknot

                # compute new fit removing deviant points (with fixed knots)
                xnewfit = []
                ynewfit = []
                for idum in range(len(xaxis1)):
                    delta_sp = abs(sp_collapsed[idum] - sp_median[idum])
                    rms_tmp = rms_array[iknot_array[idum]]
                    if idum == 0 or idum == (len(xaxis1) - 1):
                        lok = True
                    elif rms_tmp > 0:
                        if delta_sp < 3.0 * rms_tmp:
                            lok = True
                        else:
                            lok = False
                    else:
                        lok = True
                    if lok:
                        xnewfit.append(xaxis1[idum])
                        ynewfit.append(sp_collapsed[idum])
                nremove = 5
                splnew = AdaptiveLSQUnivariateSpline(
                    x=xnewfit[nremove:-nremove],
                    y=ynewfit[nremove:-nremove],
                    t=xknots[1:-1],
                    adaptive=False
                )
                sp_median = splnew(xaxis1)
            """

            ymax_spmedian = sp_median.max()
            y_threshold = ymax_spmedian * args.minimum_fraction
            lremove = np.where(sp_median < y_threshold)
            sp_median[lremove] = 0.0
            sp_mask[lremove] = True

            image2d_sp_median[islitlet - 1, :] = sp_median
            image2d_sp_mask[islitlet - 1, :] = sp_mask

            if abs(args.debugplot) % 10 != 0:
                xaxis1 = np.arange(1, naxis1_slitlet2d + 1)
                title = 'Slitlet#' + str(islitlet) + ' (median spectrum)'
                ax = ximplotxy(xaxis1, sp_collapsed,
                               title=title,
                               show=False, **{'label' : 'collapsed spectrum'})
                ax.plot(xaxis1, sp_median, label='fitted spectrum')
                ax.plot([1, naxis1_slitlet2d], 2*[y_threshold],
                        label='threshold')
                # ax.plot(xknots, yknots, 'o', label='knots')
                ax.legend()
                ax.set_ylim(-0.05*ymax_spmedian, 1.05*ymax_spmedian)
                pause_debugplot(args.debugplot,
                                pltshow=True, tight_layout=True)
        else:
            if args.debugplot == 0:
                islitlet_progress(islitlet, EMIR_NBARS, ignore=True)

    # ToDo: compute "average" spectrum for each pseudo-longslit, scaling
    #       with the median signal in each slitlet; derive a particular
    #       spectrum for each slitlet (scaling properly)

    image2d_sp_median_masked = np.ma.masked_array(
        image2d_sp_median,
        mask=image2d_sp_mask
    )
    ycut_median = np.ma.median(image2d_sp_median_masked, axis=1).data
    ycut_median_2d = np.repeat(ycut_median, EMIR_NAXIS1).reshape(
        EMIR_NBARS, EMIR_NAXIS1)
    image2d_sp_median_eq = image2d_sp_median_masked / ycut_median_2d
    image2d_sp_median_eq = image2d_sp_median_eq.data

    if True:
        ximshow(image2d_sp_median, title='sp_median', debugplot=12)
        ximplotxy(np.arange(1, EMIR_NBARS + 1), ycut_median, 'ro',
                  title='median value of each spectrum', debugplot=12)
        ximshow(image2d_sp_median_eq, title='sp_median_eq', debugplot=12)

    csu_conf_fitsfile.display_pseudo_longslits(
        list_valid_slitlets=list_valid_islitlets)
    dict_longslits = csu_conf_fitsfile.pseudo_longslits()

    # compute median spectrum for each longslit and insert (properly
    # scaled) that spectrum in each slitlet belonging to that longslit
    image2d_sp_median_longslit = np.zeros((EMIR_NBARS, EMIR_NAXIS1))
    islitlet = 1
    loop = True
    while loop:
        if islitlet in list_valid_islitlets:
            imin = dict_longslits[islitlet].imin()
            imax = dict_longslits[islitlet].imax()
            print('--> imin, imax: ', imin, imax)
            sp_median_longslit = np.median(
                image2d_sp_median_eq[(imin - 1):imax, :], axis=0)
            for i in range(imin, imax+1):
                print('----> i: ', i)
                image2d_sp_median_longslit[(i - 1), :] = \
                    sp_median_longslit * ycut_median[i - 1]
            islitlet = imax
        else:
            print('--> ignoring: ', islitlet)
        if islitlet == EMIR_NBARS:
            loop = False
        else:
            islitlet += 1
    if True:
        ximshow(image2d_sp_median_longslit, debugplot=12)

    # initialize rectified image
    image2d_flatfielded = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

    # main loop
    for islitlet in list(range(1, EMIR_NBARS + 1)):
        if islitlet in list_valid_islitlets:
            if args.debugplot == 0:
                islitlet_progress(islitlet, EMIR_NBARS, ignore=False)
            # define Slitlet2D object
            slt = Slitlet2D(islitlet=islitlet,
                            rectwv_coeff=rectwv_coeff,
                            debugplot=args.debugplot)

            # extract (distorted) slitlet from the initial image
            slitlet2d = slt.extract_slitlet2d(
                image_2k2k=image2d,
                subtitle='original image'
            )

            # rectify slitlet
            slitlet2d_rect = slt.rectify(
                slitlet2d=slitlet2d,
                resampling=args.resampling,
                subtitle='original rectified'
            )
            naxis2_slitlet2d, naxis1_slitlet2d = slitlet2d_rect.shape

            sp_median = image2d_sp_median_longslit[islitlet - 1, :]

            # generate rectified slitlet region filled with the median spectrum
            slitlet2d_rect_spmedian = np.tile(sp_median, (naxis2_slitlet2d, 1))
            if abs(args.debugplot) > 10:
                slt.ximshow_rectified(
                    slitlet2d_rect=slitlet2d_rect_spmedian,
                    subtitle='rectified, filled with median spectrum'
                )

            # unrectified image
            slitlet2d_unrect_spmedian = slt.rectify(
                slitlet2d=slitlet2d_rect_spmedian,
                resampling=args.resampling,
                inverse=True,
                subtitle='unrectified, filled with median spectrum'
            )

            # normalize initial slitlet image (avoid division by zero)
            slitlet2d_norm = np.zeros_like(slitlet2d)
            for j in range(naxis1_slitlet2d):
                for i in range(naxis2_slitlet2d):
                    den = slitlet2d_unrect_spmedian[i, j]
                    if den == 0:
                        slitlet2d_norm[i, j] = 1.0
                    else:
                        slitlet2d_norm[i, j] = slitlet2d[i, j] / den

            if abs(args.debugplot) > 10:
                slt.ximshow_unrectified(
                    slitlet2d=slitlet2d_norm,
                    subtitle='unrectified, pixel-to-pixel'
                )

            # check for pseudo-longslit with previous slitlet
            if islitlet > 1:
                if (islitlet - 1) in list_valid_islitlets:
                    c1 = csu_conf_fitsfile.csu_bar_slit_center(islitlet - 1)
                    w1 = csu_conf_fitsfile.csu_bar_slit_width(islitlet - 1)
                    c2 = csu_conf_fitsfile.csu_bar_slit_center(islitlet)
                    w2 = csu_conf_fitsfile.csu_bar_slit_width(islitlet)
                    if abs(w1-w2)/w1 < 0.25:
                        wmean = (w1 + w2) / 2.0
                        if abs(c1 - c2) < wmean/4.0:
                            same_slitlet_below = True
                        else:
                            same_slitlet_below = False
                    else:
                        same_slitlet_below = False
                else:
                    same_slitlet_below = False
            else:
                same_slitlet_below = False

            # check for pseudo-longslit with next slitlet
            if islitlet < EMIR_NBARS:
                if (islitlet + 1) in list_valid_islitlets:
                    c1 = csu_conf_fitsfile.csu_bar_slit_center(islitlet)
                    w1 = csu_conf_fitsfile.csu_bar_slit_width(islitlet)
                    c2 = csu_conf_fitsfile.csu_bar_slit_center(islitlet + 1)
                    w2 = csu_conf_fitsfile.csu_bar_slit_width(islitlet + 1)
                    if abs(w1-w2)/w1 < 0.25:
                        wmean = (w1 + w2) / 2.0
                        if abs(c1 - c2) < wmean/4.0:
                            same_slitlet_above = True
                        else:
                            same_slitlet_above = False
                    else:
                        same_slitlet_above = False
                else:
                    same_slitlet_above = False
            else:
                same_slitlet_above = False

            for j in range(EMIR_NAXIS1):
                xchannel = j + 1
                y0_lower = slt.list_frontiers[0](xchannel)
                y0_upper = slt.list_frontiers[1](xchannel)
                n1, n2 = nscan_minmax_frontiers(y0_frontier_lower=y0_lower,
                                                y0_frontier_upper=y0_upper,
                                                resize=True)
                # note that n1 and n2 are scans (ranging from 1 to NAXIS2)
                nn1 = n1 - slt.bb_ns1_orig + 1
                nn2 = n2 - slt.bb_ns1_orig + 1
                image2d_flatfielded[(n1 - 1):n2, j] = \
                    slitlet2d_norm[(nn1 - 1):nn2, j]

                # force to 1.0 region around frontiers
                if not same_slitlet_below:
                    image2d_flatfielded[(n1 - 1):(n1 + 2), j] = 1
                if not same_slitlet_above:
                    image2d_flatfielded[(n2 - 5):n2, j] = 1
        else:
            if args.debugplot == 0:
                islitlet_progress(islitlet, EMIR_NBARS, ignore=True)

    if args.debugplot == 0:
        print('OK!')

    # restore global offsets
    image2d_flatfielded = apply_integer_offsets(
        image2d=image2d_flatfielded ,
        offx=-rectwv_coeff.global_integer_offset_x_pix,
        offy=-rectwv_coeff.global_integer_offset_y_pix
    )

    # set pixels below minimum value to 1.0
    filtered = np.where(image2d_flatfielded < args.minimum_value_in_output)
    image2d_flatfielded[filtered] = 1.0

    # set pixels above maximum value to 1.0
    filtered = np.where(image2d_flatfielded > args.maximum_value_in_output)
    image2d_flatfielded[filtered] = 1.0

    # save output file
    save_ndarray_to_fits(
        array=image2d_flatfielded,
        file_name=args.outfile,
        main_header=header,
        overwrite=True
    )
    print('>>> Saving file ' + args.outfile.name)


if __name__ == "__main__":
    main()
