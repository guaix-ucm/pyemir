#
# Copyright 2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Spectroscopy mode, compute pixel-to-pixel flatfield
"""

import astropy.io.fits as fits
import contextlib
import numpy as np
import uuid

import numina.array.combine as combine
from numina.core import Parameter
from numina.core import Result
from numina.processing.combine import combine_imgs

from emirdrp.core.recipe import EmirRecipe
import emirdrp.datamodel
from emirdrp.processing.wavecal.rectwv_coeff_from_mos_library \
    import rectwv_coeff_from_mos_library
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import save_four_ds9
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 \
    import save_spectral_lines_ds9
import emirdrp.products as prods
import emirdrp.requirements as reqs

from emirdrp.core import EMIR_MINIMUM_SLITLET_WIDTH_MM
from emirdrp.core import EMIR_MAXIMUM_SLITLET_WIDTH_MM


class SpecFlatPix2Pix(EmirRecipe):
    """Process continuum exposures of continuum lamp (lamp ON-OFF)

    """

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_rectwv = reqs.MasterRectWaveRequirement(optional=True)
    rectwv_coeff = reqs.RectWaveCoeffRequirement(optional=True)

    # note that 'sum' is not allowed as combination method
    method = Parameter(
        'mean',
        description='Combination method',
        choices=['mean', 'median', 'sigmaclip']
    )
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method'
    )
    minimum_slitlet_width_mm = Parameter(
        float(EMIR_MINIMUM_SLITLET_WIDTH_MM),
        description='Minimum width (mm) for a valid slitlet',
        optional=True
    )
    maximum_slitlet_width_mm = Parameter(
        float(EMIR_MAXIMUM_SLITLET_WIDTH_MM),
        description='Maximum width (mm) for a valid slitlet',
        optional=True
    )
    global_integer_offset_x_pix = Parameter(
        0,
        description='Global offset (pixels) in wavelength direction (integer)',
        optional=True
    )
    global_integer_offset_y_pix = Parameter(
        0,
        description='Global offset (pixels) in spatial direction (integer)',
        optional=True
    )

    reduced_flatpix2pix = Result(prods.MasterSpectralFlat)

    def run(self, rinput):

        self.logger.info('starting generation of flatpix2pix')

        self.logger.info('master_rectwv.........................: {}'.format(
            rinput.master_rectwv))
        self.logger.info('rectwv_coeff..........................: {}'.format(
            rinput.rectwv_coeff))
        self.logger.info('Minimum slitlet width (mm)............: {}'.format(
            rinput.minimum_slitlet_width_mm))
        self.logger.info('Maximum slitlet width (mm)............: {}'.format(
            rinput.maximum_slitlet_width_mm))
        self.logger.info('Global offset X direction (pixels)....: {}'.format(
            rinput.global_integer_offset_x_pix))
        self.logger.info('Global offset Y direction (pixels)....: {}'.format(
            rinput.global_integer_offset_y_pix))

        # check rectification and wavelength calibration information
        if rinput.master_rectwv is None and rinput.rectwv_coeff is None:
            raise ValueError('No master_rectwv nor rectwv_coeff data have '
                             'been provided')
        elif rinput.master_rectwv is not None and \
                rinput.rectwv_coeff is not None:
            raise ValueError('master_rectwv and rectwv_coeff cannot be '
                             'use simultaneously')

        # check headers to detect lamp status (on/off)
        list_lampincd = []
        for fname in rinput.obresult.frames:
            with fname.open() as f:
                list_lampincd.append(f[0].header['lampincd'])

        # check number of images
        nimages = len(rinput.obresult.frames)
        n_on = list_lampincd.count(1)
        n_off = list_lampincd.count(0)
        self.logger.info('Number of images with lamp ON.........: {}'.format(
            n_on))
        self.logger.info('Number of images with lamp OFF........: {}'.format(
            n_off))
        self.logger.info('Total number of images................: {}'.format(
            nimages))
        if n_on == 0:
            raise ValueError('Insufficient number of images with lamp ON')
        if n_on + n_off != nimages:
            raise ValueError('Number of images does not match!')

        # check combination method
        if rinput.method != 'sigmaclip':
            if rinput.method_kwargs != {}:
                raise ValueError('Unexpected method_kwargs={}'.format(
                    rinput.method_kwargs))

        # build object to proceed with bpm, bias, and dark (not flat)
        flow = self.init_filters(rinput)

        # available combination methods
        fmethod = getattr(combine, rinput.method)

        # basic reduction of images with lamp ON or OFF
        lampmode = {0: 'off', 1: 'on'}
        reduced_image_on = None
        reduced_image_off = None
        for imode in lampmode.keys():
            self.logger.info('starting basic reduction of images with'
                             ' lamp {}'.format(lampmode[imode]))
            tmplist = [rinput.obresult.frames[i] for i, lampincd in
                       enumerate(list_lampincd) if lampincd == imode]
            if len(tmplist) > 0:
                with contextlib.ExitStack() as stack:
                    hduls = [stack.enter_context(fname.open()) for fname
                             in tmplist]
                    reduced_image = combine_imgs(
                        hduls,
                        method=fmethod,
                        method_kwargs=rinput.method_kwargs,
                        errors=False,
                        prolog=None
                    )
                if imode == 0:
                    reduced_image_off = flow(reduced_image)
                    hdr = reduced_image_off[0].header
                    self.set_base_headers(hdr)
                    self.save_intermediate_img(reduced_image_off,
                                               'reduced_image_off.fits')
                elif imode == 1:
                    reduced_image_on = flow(reduced_image)
                    hdr = reduced_image_on[0].header
                    self.set_base_headers(hdr)
                    self.save_intermediate_img(reduced_image_on,
                                               'reduced_image_on.fits')
                else:
                    raise ValueError('Unexpected imode={}'.format(imode))

        # computation of ON-OFF
        header_on = reduced_image_on[0].header
        data_on = reduced_image_on[0].data.astype('float32')
        if n_off > 0:
            header_off = reduced_image_off[0].header
            data_off = reduced_image_off[0].data.astype('float32')
        else:
            header_off = None
            data_off = np.zeros_like(data_on)
        reduced_data = data_on - data_off

        # update reduced image header
        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(fname.open()) for fname in
                     rinput.obresult.frames]
            reduced_image = self.create_reduced_image(
                hduls,
                reduced_data,
                header_on,
                header_off,
                list_lampincd,
                header_mos_onoff=None,
            )

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # define rectification and wavelength calibration coefficients
        if rinput.rectwv_coeff is None:
            rectwv_coeff = rectwv_coeff_from_mos_library(
                reduced_image,
                rinput.master_rectwv
            )
        else:
            rectwv_coeff = rinput.rectwv_coeff
        # set global offsets
        rectwv_coeff.global_integer_offset_x_pix = \
            rinput.global_integer_offset_x_pix
        rectwv_coeff.global_integer_offset_y_pix = \
            rinput.global_integer_offset_y_pix
        # save as JSON in work directory
        self.save_structured_as_json(rectwv_coeff, 'rectwv_coeff.json')
        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(rectwv_coeff)
            save_spectral_lines_ds9(rectwv_coeff)

        # ToDo: continue here
        # Include code from tools/continuum_flatfield.py

    def create_reduced_image(self, hduls, reduced_data,
                             header_on, header_off,
                             list_lampincd,
                             header_mos_onoff):
        # Copy header of first image
        base_header = hduls[0][0].header.copy()

        hdu = fits.PrimaryHDU(reduced_data, header=base_header)
        self.set_base_headers(hdu.header)

        self.logger.debug('update result header')
        if header_mos_onoff is not None:
            self.logger.debug('update result header')
            crpix1 = header_mos_onoff['crpix1']
            crval1 = header_mos_onoff['crval1']
            cdelt1 = header_mos_onoff['cdelt1']

            # update wavelength calibration in FITS header
            for keyword in ['crval1', 'crpix1', 'crval2', 'crpix2']:
                if keyword in hdu.header:
                    hdu.header.remove(keyword)
            hdu.header['crpix1'] = (crpix1, 'reference pixel')
            hdu.header['crval1'] = (crval1, 'central wavelength at crpix1')
            hdu.header['cdelt1'] = \
                (cdelt1, 'linear dispersion (Angstrom/pixel)')
            hdu.header['cunit1'] = ('Angstrom', 'units along axis1')
            hdu.header['ctype1'] = 'WAVELENGTH'
            hdu.header['crpix2'] = (0.0, 'reference pixel')
            hdu.header['crval2'] = (0.0, 'central value at crpix2')
            hdu.header['cdelt2'] = (1.0, 'increment')
            hdu.header['ctype2'] = 'PIXEL'
            hdu.header['cunit2'] = ('Pixel', 'units along axis2')
            for keyword in ['cd1_1', 'cd1_2', 'cd2_1', 'cd2_2',
                            'PCD1_1', 'PCD1_2', 'PCD2_1', 'PCD2_2',
                            'PCRPIX1', 'PCRPIX2']:
                if keyword in hdu.header:
                    hdu.header.remove(keyword)

        # update additional keywords
        hdu.header['UUID'] = str(uuid.uuid1())
        hdu.header['OBSMODE'] = 'flatpix2pix'
        hdu.header['TSUTC2'] = hduls[-1][0].header['TSUTC2']
        hdu.header['history'] = "Processed flatpix2pix"
        hdu.header['NUM-NCOM'] = (len(hduls), 'Number of combined frames')

        # update history
        dm = emirdrp.datamodel.EmirDataModel()
        for img, lampincd in zip(hduls, list_lampincd):
            imgid = dm.get_imgid(img)
            hdu.header['HISTORY'] = "Image '{}' has lampincd='{}'".format(
                imgid, lampincd)
        hdu.header['HISTORY'] = "Processed flatpix2pix"
        hdu.header['HISTORY'] = '--- Reduction of images with lamp ON ---'
        for line in header_on['HISTORY']:
            hdu.header['HISTORY'] = line
        if header_off is not None:
            hdu.header['HISTORY'] = '--- Reduction of images with lamp OFF ---'
            for line in header_off['HISTORY']:
                hdu.header['HISTORY'] = line
        result = fits.HDUList([hdu])
        return result

    def set_base_headers(self, hdr):
        newhdr = super(SpecFlatPix2Pix, self).set_base_headers(hdr)
        # Update EXP to 0
        newhdr['EXP'] = 0
        return newhdr
