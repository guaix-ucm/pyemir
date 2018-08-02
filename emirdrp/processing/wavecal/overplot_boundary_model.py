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
import os.path
import sys

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow
from numina.tools.arg_file_is_new import arg_file_is_new

from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.products import MasterRectWave
from emirdrp.tools.fit_boundaries import bound_params_from_dict
from emirdrp.tools.fit_boundaries import overplot_boundaries_from_params
from emirdrp.tools.fit_boundaries import overplot_frontiers_from_params
from emirdrp.tools.fit_boundaries import save_boundaries_from_params_ds9
from emirdrp.tools.fit_boundaries import save_frontiers_from_params_ds9

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NBARS


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
    parser.add_argument("--ds9_boundaries",
                        help="Output ds9 region file with slitlet boundaries",
                        type=lambda x: arg_file_is_new(parser, x))
    parser.add_argument("--ds9_frontiers",
                        help="Output ds9 region file with slitlet frontiers",
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

    # read input FITS file
    hdulist = fits.open(args.fitsfile)
    image_header = hdulist[0].header
    image2d = hdulist[0].data
    hdulist.close()

    naxis1 = image_header['naxis1']
    naxis2 = image_header['naxis2']

    if image2d.shape != (naxis2, naxis1):
        raise ValueError("Unexpected error with NAXIS1, NAXIS2")

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
            global_offset_x_pix = -args.global_integer_offset_x_pix,
            global_offset_y_pix = -args.global_integer_offset_y_pix
        )

        # show plot
        pause_debugplot(12, pltshow=True)


if __name__ == "__main__":

    main()
