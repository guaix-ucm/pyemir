#
# Copyright 2018-2020 Universidad Complutense de Madrid
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
import logging
import numpy as np
import sys

from numina.array.display.logging_from_debugplot import logging_from_debugplot
from numina.array.wavecalib.apply_integer_offsets import apply_integer_offsets
from numina.array.wavecalib.fix_pix_borders import find_pix_borders
from numina.array.wavecalib.resample import resample_image2d_flux
from numina.frame.utils import copy_img
from numina.tools.arg_file_is_new import arg_file_is_new

import emirdrp.datamodel as datamodel
from emirdrp.instrument.dtuconf import DtuConf
from emirdrp.products import RectWaveCoeff

from .set_wv_parameters import set_wv_parameters
from .slitlet2d import Slitlet2D

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_NPIXPERSLIT_RECTIFIED


def apply_rectwv_coeff(reduced_image,
                       rectwv_coeff,
                       args_resampling=2,
                       args_ignore_dtu_configuration=True,
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
    args_resampling : int
        1: nearest neighbour, 2: flux preserving interpolation.
    args_ignore_dtu_configuration : bool
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

    rectwv_image = copy_img(reduced_image)
    header = rectwv_image[0].header
    image2d = rectwv_image[0].data

    # apply global offsets
    image2d = apply_integer_offsets(
        image2d=image2d,
        offx=rectwv_coeff.global_integer_offset_x_pix,
        offy=rectwv_coeff.global_integer_offset_y_pix
    )

    # check grism and filter
    filter_name = header['filter']
    logger.info('Filter: ' + filter_name)
    if filter_name != rectwv_coeff.tags['filter']:
        raise ValueError('Filter name does not match!')
    grism_name = header['grism']
    logger.info('Grism: ' + grism_name)
    if grism_name != rectwv_coeff.tags['grism']:
        raise ValueError('Grism name does not match!')

    # read the DTU configuration from the image header
    mecs_header = datamodel.get_mecs_header(reduced_image)
    dtu_conf = DtuConf.from_img(reduced_image)

    # retrieve DTU configuration from RectWaveCoeff object
    dtu_conf_calib = DtuConf.from_values(
        **rectwv_coeff.meta_info['dtu_configuration']
    )
    # check that the DTU configuration employed to obtain the calibration
    # corresponds to the DTU configuration in the input FITS file
    if dtu_conf != dtu_conf_calib:
        if args_ignore_dtu_configuration:
            logger.warning('DTU configuration differences found!')
        else:
            logger.warning('DTU configuration from image header:')
            logger.warning(dtu_conf)
            logger.warning('DTU configuration from master calibration:')
            logger.warning(dtu_conf_calib)
            raise ValueError("DTU configurations do not match!")
    else:
        logger.info('DTU configuration match!')

    # valid slitlet numbers
    list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
    for idel in rectwv_coeff.missing_slitlets:
        list_valid_islitlets.remove(idel)
    logger.debug('Valid slitlet numbers:\n' + str(list_valid_islitlets))

    # ---

    # relevant wavelength calibration parameters for rectified and wavelength
    # calibrated image
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    crpix1_enlarged = wv_parameters['crpix1_enlarged']
    crval1_enlarged = wv_parameters['crval1_enlarged']
    cdelt1_enlarged = wv_parameters['cdelt1_enlarged']
    naxis1_enlarged = wv_parameters['naxis1_enlarged']

    # initialize rectified and wavelength calibrated image
    naxis2_enlarged = EMIR_NBARS * EMIR_NPIXPERSLIT_RECTIFIED
    image2d_rectwv = np.zeros((naxis2_enlarged, naxis1_enlarged),
                              dtype='float32')

    # main loop

    logger.info('Applying rectification and wavelength calibration')
    logger.info('RectWaveCoeff uuid={}'.format(rectwv_coeff.uuid))

    cout = '0'
    for islitlet in range(1, EMIR_NBARS + 1):

        if islitlet in list_valid_islitlets:

            # define Slitlet2D object
            slt = Slitlet2D(islitlet=islitlet,
                            rectwv_coeff=rectwv_coeff,
                            debugplot=debugplot)

            # extract (distorted) slitlet from the initial image
            slitlet2d = slt.extract_slitlet2d(image2d)

            # rectify slitlet
            slitlet2d_rect = slt.rectify(slitlet2d, resampling=args_resampling)

            # wavelength calibration of the rectifed slitlet
            slitlet2d_rect_wv = resample_image2d_flux(
                image2d_orig=slitlet2d_rect,
                naxis1=naxis1_enlarged,
                cdelt1=cdelt1_enlarged,
                crval1=crval1_enlarged,
                crpix1=crpix1_enlarged,
                coeff=slt.wpoly
            )

            # minimum and maximum useful row in the full 2d rectified image
            # (starting from 0)
            i1 = slt.iminslt - 1
            i2 = slt.imaxslt

            # minimum and maximum scan in the rectified slitlet
            # (in pixels, from 1 to NAXIS2)
            ii1 = slt.min_row_rectified
            ii2 = slt.max_row_rectified + 1

            # save rectified slitlet in its corresponding location within
            # the full 2d rectified image
            image2d_rectwv[i1:i2, :] = slitlet2d_rect_wv[ii1:ii2, :]

            # include scan range in FITS header
            header['imnslt' + str(islitlet).zfill(2)] = \
                slt.iminslt, 'minimum Y pixel of useful slitlet region'
            header['imxslt' + str(islitlet).zfill(2)] = \
                slt.imaxslt, 'maximum Y pixel of useful slitlet region'

            # determine useful channel region in each spectrum and include
            # that information in FITS header
            jminslt = []
            jmaxslt = []
            for idum in range(ii1, ii2 + 1):
                jminmax = find_pix_borders(
                    slitlet2d_rect_wv[idum, :],
                    sought_value=0
                )
                if jminmax != (-1, naxis1_enlarged):
                    jminslt.append(jminmax[0])
                    jmaxslt.append(jminmax[1])
            if len(jminslt) > 0:
                slt.jminslt = min(jminslt) + 1
                slt.jmaxslt = max(jmaxslt) + 1
            header['jmnslt' + str(islitlet).zfill(2)] = \
                slt.jminslt, 'minimum X pixel of useful slitlet region'
            header['jmxslt' + str(islitlet).zfill(2)] = \
                slt.jmaxslt, 'maximum X pixel of useful slitlet region'

            cout += '.'

        else:

            # include scan and channel range in FITS header
            header['imnslt' + str(islitlet).zfill(2)] = \
                0, 'minimum Y pixel of useful slitlet region'
            header['imxslt' + str(islitlet).zfill(2)] = \
                0, 'maximum Y pixel of useful slitlet region'
            header['jmnslt' + str(islitlet).zfill(2)] = \
                0, 'minimum X pixel of useful slitlet region'
            header['jmxslt' + str(islitlet).zfill(2)] = \
                0, 'maximum X pixel of useful slitlet region'

            cout += 'i'

        if islitlet % 10 == 0:
            if cout != 'i':
                cout = str(islitlet // 10)

        logger.info(cout)

    # update wavelength calibration in FITS header
    logger.info('Updating image header')
    for keyword in ['crval1', 'crpix1', 'crval2', 'crpix2']:
        if keyword in header:
            header.remove(keyword)
    header['crpix1'] = (crpix1_enlarged, 'reference pixel')
    header['crval1'] = (crval1_enlarged, 'central wavelength at crpix1')
    header['cdelt1'] = (cdelt1_enlarged, 'linear dispersion (Angstrom/pixel)')
    header['cunit1'] = ('Angstrom', 'units along axis1')
    header['ctype1'] = 'WAVE'
    header['crpix2'] = (0.0, 'reference pixel')
    header['crval2'] = (0.0, 'central value at crpix2')
    header['cdelt2'] = (1.0, 'increment')
    header['ctype2'] = 'PIXEL'
    header['cunit2'] = ('Pixel', 'units along axis2')
    for keyword in ['cd1_1', 'cd1_2', 'cd2_1', 'cd2_2',
                    'PCD1_1', 'PCD1_2', 'PCD2_1', 'PCD2_2',
                    'PCRPIX1', 'PCRPIX2']:
        if keyword in header:
            header.remove(keyword)

    # update history in FITS header
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
    rectwv_image[0].data = image2d_rectwv
    return rectwv_image

    rectwv_image[0].data = image2d_rectwv
    rectwv_image[0].header = header
    return rectwv_image


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: apply rectification and wavelength '
                    'calibration polynomials for the CSU configuration of a '
                    'particular image'
    )

    # required arguments
    parser.add_argument("fitsfile",
                        help="Input FITS file",
                        type=argparse.FileType('rb'))
    parser.add_argument("--rectwv_coeff", required=True,
                        help="Input JSON file with rectification and "
                             "wavelength calibration coefficients",
                        type=argparse.FileType('rt'))
    parser.add_argument("--outfile", required=True,
                        help="Output FITS file with rectified and "
                             "wavelength calibrated image",
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
    parser.add_argument("--ignore_dtu_configuration",
                        help="Ignore DTU configurations differences between "
                             "transformation and input image",
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

    # ---

    logging_from_debugplot(args.debugplot)

    # generate RectWaveCoeff object
    rectwv_coeff = RectWaveCoeff._datatype_load(args.rectwv_coeff.name)

    # modify (when requested) global offsets
    rectwv_coeff.global_integer_offset_x_pix += \
        args.delta_global_integer_offset_x_pix
    rectwv_coeff.global_integer_offset_y_pix += \
        args.delta_global_integer_offset_y_pix

    # generate HDUList object
    # read FITS image and its corresponding header
    hdulist = fits.open(args.fitsfile)

    # rectification and wavelength calibration
    reduced_arc = apply_rectwv_coeff(
        hdulist,
        rectwv_coeff,
        args_resampling=args.resampling,
        args_ignore_dtu_configuration=args.ignore_dtu_configuration,
        debugplot=args.debugplot
    )

    # save result
    reduced_arc.writeto(args.outfile, overwrite=True)


if __name__ == "__main__":
    main()
