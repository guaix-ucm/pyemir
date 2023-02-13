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
import os.path
import pkgutil
from six import StringIO
import sys

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow
from numina.array.distortion import fmap
from numina.tools.arg_file_is_new import arg_file_is_new

from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.products import MasterRectWave
from emirdrp.processing.wavecal.rectwv_coeff_from_mos_library \
    import rectwv_coeff_from_mos_library
from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.fit_boundaries import overplot_boundaries_from_params
from emirdrp.tools.fit_boundaries import overplot_frontiers_from_params
from emirdrp.tools.fit_boundaries import save_boundaries_from_params_ds9
from emirdrp.tools.fit_boundaries import save_frontiers_from_params_ds9

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS
from emirdrp.tools.fit_boundaries import EXPECTED_PARAMETER_LIST
from emirdrp.tools.fit_boundaries import EXPECTED_PARAMETER_LIST_EXTENDED


def overplot_lines(ax,
                   catlines_all_wave,
                   list_valid_islitlets,
                   rectwv_coeff,
                   global_integer_offset_x_pix, global_integer_offset_y_pix,
                   ds9_file, debugplot):
    """Overplot lines (arc/OH).

    Parameters
    ----------
    ax : matplotlib axes
        Current plot axes.
    catlines_all_wave : numpy array
        Array with wavelengths of the lines to be overplotted.
    list_valid_islitlets : list of integers
        List with numbers of valid slitlets.
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    global_integer_offset_x_pix : int
        Global offset in the X direction to be applied after computing
        the expected location.
    global_integer_offset_y_pix : int
        Global offset in the Y direction to be applied after computing
        the expected location.
    ds9_file : file handler or None
        File handler to ds9 region file where the location of lines
        must be saved.
    debugplot : int
        Debugging level for messages and plots. For details see
        'numina.array.display.pause_debugplot.py'.

    """

    for islitlet in list_valid_islitlets:
        crval1_linear = rectwv_coeff.contents[islitlet - 1]['crval1_linear']
        cdelt1_linear = rectwv_coeff.contents[islitlet - 1]['cdelt1_linear']
        crvaln_linear = crval1_linear + (EMIR_NAXIS1 - 1) * cdelt1_linear
        bb_ns1_orig = rectwv_coeff.contents[islitlet - 1]['bb_ns1_orig']
        ttd_order = rectwv_coeff.contents[islitlet - 1]['ttd_order']
        aij = rectwv_coeff.contents[islitlet - 1]['ttd_aij']
        bij = rectwv_coeff.contents[islitlet - 1]['ttd_bij']
        min_row_rectified = float(
            rectwv_coeff.contents[islitlet - 1]['min_row_rectified']
        )
        max_row_rectified = float(
            rectwv_coeff.contents[islitlet - 1]['max_row_rectified']
        )
        mean_row_rectified = (min_row_rectified + max_row_rectified) / 2
        wpoly_coeff = rectwv_coeff.contents[islitlet - 1]['wpoly_coeff']
        x0 = []
        y0 = []
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        for line in catlines_all_wave:
            if crval1_linear <= line <= crvaln_linear:
                tmp_coeff = np.copy(wpoly_coeff)
                tmp_coeff[0] -= line
                tmp_xroots = np.polynomial.Polynomial(tmp_coeff).roots()
                for dum in tmp_xroots:
                    if np.isreal(dum):
                        dum = dum.real
                        if 1 <= dum <= EMIR_NAXIS1:
                            x0.append(dum)
                            y0.append(mean_row_rectified)
                            x1.append(dum)
                            y1.append(min_row_rectified)
                            x2.append(dum)
                            y2.append(max_row_rectified)
        xx0, yy0 = fmap(ttd_order, aij, bij, np.array(x0), np.array(y0))
        xx0 -= global_integer_offset_x_pix
        yy0 += bb_ns1_orig
        yy0 -= global_integer_offset_y_pix
        xx1, yy1 = fmap(ttd_order, aij, bij, np.array(x1), np.array(y1))
        xx1 -= global_integer_offset_x_pix
        yy1 += bb_ns1_orig
        yy1 -= global_integer_offset_y_pix
        xx2, yy2 = fmap(ttd_order, aij, bij, np.array(x2), np.array(y2))
        xx2 -= global_integer_offset_x_pix
        yy2 += bb_ns1_orig
        yy2 -= global_integer_offset_y_pix

        if abs(debugplot) % 10 != 0:
            if abs(debugplot) == 22:
                for xx1_, xx2_, yy1_, yy2_ in zip(xx1, xx2, yy1, yy2):
                    ax.plot([xx1_, xx2_], [yy1_, yy2_], 'c-', linewidth=2.0)
            else:
                ax.plot(xx0, yy0, 'c.')

        if ds9_file is not None:
            ds9_file.write(
                '#\n# islitlet...........: {0}\n'.format(islitlet)
            )
            for xx0_, yy0_ in zip(xx0, yy0):
                ds9_file.write(
                    'circle {0} {1} 2 # fill=1\n'.format(
                        xx0_, yy0_)
                )


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: overplot boundary model over FITS image'
    )

    # positional arguments
    parser.add_argument("fitsfile",
                        help="FITS file name to be displayed",
                        type=argparse.FileType('rb'))
    parser.add_argument("--rect_wpoly_MOSlibrary", required=True,
                        help="Input JSON file with library of rectification "
                             "and wavelength calibration coefficients",
                        type=argparse.FileType('rt'))

    # optional arguments
    parser.add_argument("--global_integer_offset_x_pix",
                        help="Global integer offset in the X direction "
                             "(default=0)",
                        default=0, type=int)
    parser.add_argument("--global_integer_offset_y_pix",
                        help="Global integer offset in the Y direction "
                             "(default=0)",
                        default=0, type=int)
    parser.add_argument("--arc_lines",
                        help="Overplot arc lines",
                        action="store_true")
    parser.add_argument("--oh_lines",
                        help="Overplot OH lines",
                        action="store_true")
    parser.add_argument("--ds9_frontiers",
                        help="Output ds9 region file with slitlet frontiers",
                        type=lambda x: arg_file_is_new(parser, x))
    parser.add_argument("--ds9_boundaries",
                        help="Output ds9 region file with slitlet boundaries",
                        type=lambda x: arg_file_is_new(parser, x))
    parser.add_argument("--ds9_lines",
                        help="Output ds9 region file with arc/oh lines",
                        type=lambda x: arg_file_is_new(parser, x))
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting/debugging" +
                             " (default=12)",
                        type=int, default=12,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args()

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # ---

    # avoid incompatible options
    if args.arc_lines and args.oh_lines:
        raise ValueError("--arc_lines and --oh_lines cannot be used "
                         "simultaneously")

    # --ds9_lines requires --arc_lines or --oh_lines
    if args.ds9_lines:
        if not (args.arc_lines or args.oh_lines):
            raise ValueError("--ds9_lines requires the use of either "
                             "--arc_lines or --oh_lines")

    # read input FITS file
    hdulist = fits.open(args.fitsfile)
    for item in hdulist:
        print(item)
    image_header = hdulist[0].header
    image2d = hdulist[0].data
    # hdulist.close()

    naxis1 = image_header['naxis1']
    naxis2 = image_header['naxis2']

    if image2d.shape != (naxis2, naxis1):
        raise ValueError("Unexpected error with NAXIS1, NAXIS2")
    if image2d.shape != (EMIR_NAXIS2, EMIR_NAXIS1):
        raise ValueError("Unexpected values for NAXIS1, NAXIS2")

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

    # ---

    # generate MasterRectWave object
    master_rectwv = MasterRectWave._datatype_load(
        args.rect_wpoly_MOSlibrary.name
    )

    # check that grism and filter are the expected ones
    grism_ = master_rectwv.tags['grism']
    if grism_ != grism:
        raise ValueError('Unexpected grism: ' + str(grism_))
    spfilter_ = master_rectwv.tags['filter']
    if spfilter_ != spfilter:
        raise ValueError('Unexpected filter ' + str(spfilter_))

    # valid slitlet numbers
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in master_rectwv.missing_slitlets:
        list_valid_islitlets.remove(idel)

    # read CsuConfiguration object from FITS file
    csu_config = CsuConfiguration.define_from_fits(args.fitsfile)

    # list with csu_bar_slit_center for valid slitlets
    list_csu_bar_slit_center = []
    for islitlet in list_valid_islitlets:
        list_csu_bar_slit_center.append(
            csu_config.csu_bar_slit_center(islitlet)
        )

    # define parmodel and params
    fitted_bound_param_json = {
        'contents': master_rectwv.meta_info['refined_boundary_model']
    }
    parmodel = fitted_bound_param_json['contents']['parmodel']
    fitted_bound_param_json.update({'meta_info': {'parmodel': parmodel}})
    params = bound_params_from_dict(fitted_bound_param_json)
    if parmodel != "multislit":
        raise ValueError('parmodel = "multislit" not found')

    # ---

    # define lines to be overplotted
    if args.arc_lines or args.oh_lines:

        rectwv_coeff = rectwv_coeff_from_mos_library(hdulist, master_rectwv)
        rectwv_coeff.global_integer_offset_x_pix = \
            args.global_integer_offset_x_pix
        rectwv_coeff.global_integer_offset_y_pix = \
            args.global_integer_offset_y_pix
        # rectwv_coeff.writeto('xxx.json')

        if args.arc_lines:
            if grism == 'LR':
                catlines_file = 'lines_argon_neon_xenon_empirical_LR.dat'
            else:
                catlines_file = 'lines_argon_neon_xenon_empirical.dat'
            dumdata = pkgutil.get_data('emirdrp.instrument.configs',
                                       catlines_file)
            arc_lines_tmpfile = StringIO(dumdata.decode('utf8'))
            catlines = np.genfromtxt(arc_lines_tmpfile)
            # define wavelength and flux as separate arrays
            catlines_all_wave = catlines[:, 0]
            catlines_all_flux = catlines[:, 1]
        elif args.oh_lines:
            dumdata = pkgutil.get_data('emirdrp.instrument.configs',
                                       'Oliva_etal_2013.dat')
            oh_lines_tmpfile = StringIO(dumdata.decode('utf8'))
            catlines = np.genfromtxt(oh_lines_tmpfile)
            # define wavelength and flux as separate arrays
            catlines_all_wave = np.concatenate((catlines[:, 1],
                                                catlines[:, 0]))
            catlines_all_flux = np.concatenate((catlines[:, 2],
                                                catlines[:, 2]))
        else:
            raise ValueError("This should not happen!")

    else:
            rectwv_coeff = None
            catlines_all_wave = None
            catlines_all_flux = None

    # ---

    # generate output ds9 region file with slitlet boundaries
    if args.ds9_boundaries is not None:
        save_boundaries_from_params_ds9(
            params=params,
            parmodel=parmodel,
            list_islitlet=list_valid_islitlets,
            list_csu_bar_slit_center=list_csu_bar_slit_center,
            uuid=master_rectwv.uuid,
            grism=grism,
            spfilter=spfilter,
            ds9_filename=args.ds9_boundaries.name,
            global_offset_x_pix=-args.global_integer_offset_x_pix,
            global_offset_y_pix=-args.global_integer_offset_y_pix
        )

    # generate output ds9 region file with slitlet frontiers
    if args.ds9_frontiers is not None:
        save_frontiers_from_params_ds9(
            params=params,
            parmodel=parmodel,
            list_islitlet=list_valid_islitlets,
            list_csu_bar_slit_center=list_csu_bar_slit_center,
            uuid=master_rectwv.uuid,
            grism=grism,
            spfilter=spfilter,
            ds9_filename=args.ds9_frontiers.name,
            global_offset_x_pix=-args.global_integer_offset_x_pix,
            global_offset_y_pix=-args.global_integer_offset_y_pix
        )

    # ---

    # display full image
    if abs(args.debugplot) % 10 != 0:
        ax = ximshow(image2d=image2d,
                     title=sfitsfile + "\ngrism=" + grism +
                           ", filter=" + spfilter +
                           ", rotang=" + str(round(rotang, 2)),
                     image_bbox=(1, naxis1, 1, naxis2), show=False)

        # overplot boundaries
        overplot_boundaries_from_params(
            ax=ax,
            params=params,
            parmodel=parmodel,
            list_islitlet=list_valid_islitlets,
            list_csu_bar_slit_center=list_csu_bar_slit_center,
            global_offset_x_pix=-args.global_integer_offset_x_pix,
            global_offset_y_pix=-args.global_integer_offset_y_pix
        )

        # overplot frontiers
        overplot_frontiers_from_params(
            ax=ax,
            params=params,
            parmodel=parmodel,
            list_islitlet=list_valid_islitlets,
            list_csu_bar_slit_center=list_csu_bar_slit_center,
            micolors=('b', 'b'), linetype='-',
            labels=False,    # already displayed with the boundaries
            global_offset_x_pix=-args.global_integer_offset_x_pix,
            global_offset_y_pix=-args.global_integer_offset_y_pix
        )

    else:
        ax = None

    # overplot lines
    if catlines_all_wave is not None:

        if args.ds9_lines is None:
            ds9_file = None
        else:
            ds9_file = open(args.ds9_lines.name, 'w')
            ds9_file.write('# Region file format: DS9 version 4.1\n')
            ds9_file.write('global color=#00ffff dashlist=0 0 width=2 '
                           'font="helvetica 10 normal roman" select=1 '
                           'highlite=1 dash=0 fixed=0 edit=1 '
                           'move=1 delete=1 include=1 source=1\n')
            ds9_file.write('physical\n#\n')

            ds9_file.write('#\n# uuid..: {0}\n'.format(master_rectwv.uuid))
            ds9_file.write('# filter: {0}\n'.format(spfilter))
            ds9_file.write('# grism.: {0}\n'.format(grism))
            ds9_file.write('#\n# global_offset_x_pix: {0}\n'.format(
                args.global_integer_offset_x_pix))
            ds9_file.write('# global_offset_y_pix: {0}\n#\n'.format(
                args.global_integer_offset_y_pix))
            if parmodel == "longslit":
                for dumpar in EXPECTED_PARAMETER_LIST:
                    parvalue = params[dumpar].value
                    ds9_file.write('# {0}: {1}\n'.format(dumpar, parvalue))
            else:
                for dumpar in EXPECTED_PARAMETER_LIST_EXTENDED:
                    parvalue = params[dumpar].value
                    ds9_file.write('# {0}: {1}\n'.format(dumpar, parvalue))

        overplot_lines(ax,
                       catlines_all_wave,
                       list_valid_islitlets,
                       rectwv_coeff,
                       args.global_integer_offset_x_pix,
                       args.global_integer_offset_y_pix,
                       ds9_file,
                       args.debugplot
                       )

        if ds9_file is not None:
            ds9_file.close()

    if ax is not None:
        # show plot
        pause_debugplot(12, pltshow=True)


if __name__ == "__main__":

    main()
