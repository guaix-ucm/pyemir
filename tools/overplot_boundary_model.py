from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
import json
import os.path
import sys

from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.display.ximshow import ximshow
from emirdrp.instrument.csu_configuration import CsuConfiguration

from arg_file_is_new import arg_file_is_new
from fit_boundaries import bound_params_from_dict
from fit_boundaries import overplot_boundaries_from_params
from fit_boundaries import overplot_frontiers_from_params
from fit_boundaries import save_boundaries_from_params_ds9
from fit_boundaries import save_frontiers_from_params_ds9

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("fitsfile",
                        help="FITS file name to be displayed",
                        type=argparse.FileType('r'))
    parser.add_argument("--fitted_bound_param", required=True,
                        help="JSON file with fitted boundary coefficients "
                             "corresponding to the multislit model",
                        type=argparse.FileType('r'))

    # optional arguments
    parser.add_argument("--ds9reg_boundaries",
                        help="Output ds9 region file with slitlet boundaries",
                        type=lambda x: arg_file_is_new(parser, x))
    parser.add_argument("--ds9reg_frontiers",
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

    # read input FITS file
    hdulist = fits.open(args.fitsfile.name)
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

    # read fitted_bound_param JSON file
    fittedpar_dict = json.loads(open(args.fitted_bound_param.name).read())
    params = bound_params_from_dict(fittedpar_dict)
    if abs(args.debugplot) in [21, 22]:
        params.pretty_print()

    parmodel = fittedpar_dict['meta-info']['parmodel']
    if parmodel != 'multislit':
        raise ValueError("Unexpected parameter model: ", parmodel)

    # define slitlet range
    islitlet_min = fittedpar_dict['tags']['islitlet_min']
    islitlet_max = fittedpar_dict['tags']['islitlet_max']
    list_islitlet = range(islitlet_min, islitlet_max + 1)

    # read CsuConfiguration object from FITS file
    csu_config = CsuConfiguration.define_from_fits(args.fitsfile)

    # define csu_bar_slit_center associated to each slitlet
    list_csu_bar_slit_center = []
    for islitlet in list_islitlet:
        list_csu_bar_slit_center.append(
            csu_config.csu_bar_slit_center(islitlet))

    # generate output ds9 region file with slitlet boundaries
    if args.ds9reg_boundaries is not None:
        save_boundaries_from_params_ds9(
            params=params,
            parmodel=parmodel,
            list_islitlet=list_islitlet,
            list_csu_bar_slit_center=list_csu_bar_slit_center,
            uuid=fittedpar_dict['uuid'],
            grism=fittedpar_dict['tags']['grism'],
            spfilter=fittedpar_dict['tags']['filter'],
            ds9_filename=args.ds9reg_boundaries.name
        )

    # generate output ds9 region file with slitlet frontiers
    if args.ds9reg_frontiers is not None:
        save_frontiers_from_params_ds9(
            params=params,
            parmodel=parmodel,
            list_islitlet=list_islitlet,
            list_csu_bar_slit_center=list_csu_bar_slit_center,
            uuid=fittedpar_dict['uuid'],
            grism=fittedpar_dict['tags']['grism'],
            spfilter=fittedpar_dict['tags']['filter'],
            ds9_filename=args.ds9reg_frontiers.name
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
            list_islitlet=list_islitlet,
            list_csu_bar_slit_center=list_csu_bar_slit_center
        )

        # overplot frontiers
        overplot_frontiers_from_params(
            ax=ax,
            params=params,
            parmodel=parmodel,
            list_islitlet=list_islitlet,
            list_csu_bar_slit_center=list_csu_bar_slit_center,
            micolors=('b', 'b'), linetype='-',
            labels=False    # already displayed with the boundaries
        )

    # show plot
    pause_debugplot(12, pltshow=True)


if __name__ == "__main__":

    main()
