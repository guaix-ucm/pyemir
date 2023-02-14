#
# Copyright 2016-2022 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Spectroscopy mode, Stare Spectra
"""

from astropy.io import fits
from datetime import datetime
import logging
import numpy as np
from skimage.registration import phase_cross_correlation

from numina.core import Result
from numina.core import Requirement, Parameter
from numina.array.combine import sigmaclip
import numina.processing as proc
from numina.processing.combine import basic_processing_with_combination

import emirdrp.datamodel as datamodel
from emirdrp.core.recipe import EmirRecipe
from emirdrp.instrument.csu_configuration import CsuConfiguration
import emirdrp.products as prods
from emirdrp.products import RectWaveCoeff
from emirdrp.processing.wavecal.apply_rectwv_coeff import apply_rectwv_coeff
from emirdrp.processing.wavecal.median_slitlets_rectified \
    import median_slitlets_rectified
from emirdrp.processing.wavecal.rectwv_coeff_from_mos_library \
    import rectwv_coeff_from_mos_library
from emirdrp.processing.wavecal.rescale_array_z1z2 import rescale_array_to_z1z2
from emirdrp.processing.wavecal.retrieve_catlines import retrieve_catlines
from emirdrp.processing.wavecal.synthetic_lines_rawdata import \
    synthetic_lines_rawdata
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import save_four_ds9
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 \
    import save_spectral_lines_ds9
from emirdrp.processing.wavecal.refine_rectwv_coeff import refine_rectwv_coeff
import emirdrp.requirements as reqs

from emirdrp.core import EMIR_MINIMUM_SLITLET_WIDTH_MM
from emirdrp.core import EMIR_MAXIMUM_SLITLET_WIDTH_MM
from emirdrp.core import EMIR_NBARS


class StareSpectraRecipe(EmirRecipe):
    """Process images in Stare spectra mode"""

    logger = logging.getLogger(__name__)

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    master_sky = reqs.SpectralSkyRequirement(optional=True)

    stare = Result(prods.ProcessedMOS)

    def run(self, rinput):
        self.logger.info('starting stare spectra reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=sigmaclip)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        # Update EXP to 0
        hdr['EXP'] = 0

        self.logger.info('end stare spectra reduction')
        result = self.create_result(stare=hdulist)
        return result


class StareSpectraWaveRecipe(EmirRecipe):
    """Process images in Stare spectra at the GTC.

    This recipe is intended to be used at GTC. The rectification
    and wavelength calibration can computed from a model if this
    model (master_rectwv) is provided as input.

    """

    logger = logging.getLogger(__name__)

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    master_rectwv = reqs.MasterRectWaveRequirement(optional=True)
    master_sky = reqs.SpectralSkyRequirement(optional=True)

    reduced_image = Result(prods.ProcessedImage)
    reduced_mos = Result(prods.ProcessedMOS)

    def run(self, rinput):
        self.logger.info('starting reduction of stare spectra')

        self.logger.info(rinput.master_rectwv)

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # apply bpm, bias, dark and flat
        reduced_image = basic_processing_with_combination(rinput, flow,
                                                          method=sigmaclip)
        # update header with additional info
        hdr = reduced_image[0].header
        self.set_base_headers(hdr)

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # rectification and wavelength calibration (if a model has
        # been provided)
        if rinput.master_rectwv:
            # RectWaveCoeff object with rectification and wavelength
            # calibration coefficients for the particular CSU configuration
            rectwv_coeff = rectwv_coeff_from_mos_library(
                reduced_image,
                rinput.master_rectwv
            )

            # apply rectification and wavelength calibration
            stare_image = apply_rectwv_coeff(
                reduced_image,
                rectwv_coeff
            )

            # save as JSON file in work directory
            self.save_structured_as_json(rectwv_coeff, 'rectwv_coeff.json')

            # ds9 region files (to be saved in the work directory)
            if self.intermediate_results:
                save_four_ds9(rectwv_coeff)
                save_spectral_lines_ds9(rectwv_coeff)

            # compute median spectra employing the useful region of the
            # rectified image
            if self.intermediate_results:
                for imode, outfile in enumerate(['median_spectra_full',
                                                 'median_spectra_slitlets',
                                                 'median_spectrum_slitlets']):
                    median_image = median_slitlets_rectified(
                        stare_image, mode=imode
                    )
                    self.save_intermediate_img(median_image, outfile + '.fits')

            # image_wl_calibrated = True

        else:

            stare_image = reduced_image

            self.logger.info('No wavelength calibration provided')
            grism_value = hdr.get('GRISM', 'unknown')
            self.logger.debug('GRISM is %s', grism_value)
            if grism_value.lower() == 'open':
                self.logger.debug('GRISM is %s, so this seems OK', grism_value)

            # image_wl_calibrated = False

        if rinput.master_sky:
            # Sky subtraction after rectification
            msky = rinput.master_sky.open()
            # Check if images have the same size.
            # if so, go ahead
            if msky[0].data.shape != stare_image[0].data.shape:
                self.logger.warning(
                    "sky and current image don't have the same shape")
            else:
                sky_corrector = proc.SkyCorrector(
                    msky[0].data,
                    datamodel=self.datamodel,
                    calibid=self.datamodel.get_imgid(msky)
                )

                stare_image = sky_corrector(stare_image)
        else:
            self.logger.info('No sky image provided')

        # save results in results directory
        self.logger.info('end reduction of stare spectra')
        result = self.create_result(reduced_image=reduced_image,
                                    reduced_mos=stare_image)
        return result

    def set_base_headers(self, hdr):
        newhdr = super(StareSpectraWaveRecipe, self).set_base_headers(hdr)
        # Update EXP to 0
        newhdr['EXP'] = 0
        return newhdr


class GenerateRectwvCoeff(EmirRecipe):
    """Process images in Stare spectra.

    This recipe generates a rectified and wavelength calibrated
    image after applying a model (master_rectwv). This calibration
    can be refined (using refine_wavecalib_mode != 0).

    """

    logger = logging.getLogger(__name__)

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    master_rectwv = reqs.MasterRectWaveRequirement()
    refine_wavecalib_mode = Parameter(
        0,
        description='Apply wavelength calibration refinement',
        optional=True,
        choices=[0, 1, 2, 11, 12]
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
    global_integer_offsets_mode = Parameter(
        'fixed',
        description='Global integer offsets computation',
        choices=['auto', 'fixed']
    )
    global_integer_offset_x_pix = Parameter(
        0,
        description='Global integer offset (pixels) in wavelength direction',
        optional=True
    )
    global_integer_offset_y_pix = Parameter(
        0,
        description='Global integer offset (pixels) in spatial direction',
        optional=True
    )

    reduced_mos = Result(prods.ProcessedMOS)
    rectwv_coeff = Result(RectWaveCoeff)

    def run(self, rinput):
        self.logger.info('starting rect.+wavecal. reduction of stare spectra')

        self.logger.info(rinput.master_rectwv)
        self.logger.info(
            'Wavelength calibration refinement mode....: {}'.format(
                rinput.refine_wavecalib_mode))
        self.logger.info(
            'Minimum slitlet width (mm)................: {}'.format(
                rinput.minimum_slitlet_width_mm))
        self.logger.info(
            'Maximum slitlet width (mm)................: {}'.format(
                rinput.maximum_slitlet_width_mm))
        self.logger.info(
            'Global integer offsets mode...............: {}'.format(
                rinput.global_integer_offsets_mode
            )
        )
        self.logger.info(
            'Global integer offset X direction (pixels): {}'.format(
                rinput.global_integer_offset_x_pix))
        self.logger.info(
            'Global integer offset Y direction (pixels): {}'.format(
                rinput.global_integer_offset_y_pix))

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # apply bpm, bias, dark and flat
        reduced_image = basic_processing_with_combination(rinput, flow,
                                                          method=sigmaclip)
        # update header with additional info
        hdr = reduced_image[0].header
        self.set_base_headers(hdr)

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # RectWaveCoeff object with rectification and wavelength
        # calibration coefficients for the particular CSU configuration
        rectwv_coeff = rectwv_coeff_from_mos_library(
            reduced_image,
            rinput.master_rectwv
        )

        # wavelength calibration refinement
        # 0 -> no refinement
        # 1 -> apply global offset to all the slitlets (using ARC lines)
        # 2 -> apply individual offset to each slitlet (using ARC lines)
        # 11 -> apply global offset to all the slitlets (using OH lines)
        # 12 -> apply individual offset to each slitlet (using OH lines)
        if rinput.refine_wavecalib_mode != 0:
            main_header = reduced_image[0].header
            mecs_header = datamodel.get_mecs_header(reduced_image)
            # determine useful slitlets
            csu_config = CsuConfiguration.define_from_header(mecs_header)
            # segregate slitlets
            list_useful_slitlets = csu_config.widths_in_range_mm(
                minwidth=rinput.minimum_slitlet_width_mm,
                maxwidth=rinput.maximum_slitlet_width_mm
            )
            # remove missing slitlets
            if len(rectwv_coeff.missing_slitlets) > 0:
                for iremove in rectwv_coeff.missing_slitlets:
                    if iremove in list_useful_slitlets:
                        list_useful_slitlets.remove(iremove)

            list_not_useful_slitlets = [i for i in
                                        list(range(1, EMIR_NBARS + 1))
                                        if i not in list_useful_slitlets]
            self.logger.info('list of useful slitlets: {}'.format(
                list_useful_slitlets))
            self.logger.info('list of unusable slitlets: {}'.format(
                list_not_useful_slitlets))

            # retrieve arc/OH lines
            catlines_all_wave, catlines_all_flux = retrieve_catlines(
                rinput.refine_wavecalib_mode,
                main_header['grism']
            )

            # global integer offsets
            if rinput.global_integer_offsets_mode == 'auto':
                if (rinput.global_integer_offset_x_pix != 0) or \
                        (rinput.global_integer_offset_y_pix != 0):
                    raise ValueError('Global integer offsets must be zero when'
                                     ' mode=auto')

                # ToDo: include additional airglow emission lines

                self.logger.info('computing synthetic image')
                # generate synthetic image
                synthetic_raw_data = synthetic_lines_rawdata(
                    catlines_all_wave,
                    catlines_all_flux,
                    list_useful_slitlets,
                    rectwv_coeff
                )
                synthetic_raw_header = main_header.copy()
                synthetic_raw_header['DATE-OBS'] = \
                    datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
                chistory = 'Synthetic image'
                synthetic_raw_header.add_history(chistory)
                hdu = fits.PrimaryHDU(synthetic_raw_data.astype('float32'),
                                      header=synthetic_raw_header)
                synthetic_raw_image = fits.HDUList([hdu])
                if self.intermediate_results:
                    self.save_intermediate_img(synthetic_raw_image,
                                               'synthetic_raw_image.fits')

                # cross-correlation to determine global integer offsets
                # (rescaling data arrays to [0, 1] before using skimage
                # function)
                data1_rs, coef1_rs = rescale_array_to_z1z2(
                    reduced_image[0].data, (0, 1)
                )
                data2_rs, coef2_rs = rescale_array_to_z1z2(
                    synthetic_raw_data, (0, 1)
                )
                reference_mask = np.ones(data1_rs.shape, dtype=bool)
                moving_mask = (data2_rs > 0)
                # the number of arguments returned by phase_cross_correlation
                # is only 1 when reference_mask and moving_mask are used
                shifts = phase_cross_correlation(
                    data1_rs, data2_rs, upsample_factor=100,
                    reference_mask=reference_mask,
                    moving_mask=moving_mask, overlap_ratio=0.9)
                self.logger.info('global_float_offset_x_pix..: {}'.format(
                    -shifts[1]
                ))
                self.logger.info('global_float_offset_y_pix..: {}'.format(
                    -shifts[0]
                ))
                rectwv_coeff.global_integer_offset_x_pix = \
                    -int(round(shifts[1]))
                rectwv_coeff.global_integer_offset_y_pix = \
                    -int(round(shifts[0]))
                self.logger.info('global_integer_offset_x_pix: {}'.format(
                    rectwv_coeff.global_integer_offset_x_pix
                ))
                self.logger.info('global_integer_offset_y_pix: {}'.format(
                    rectwv_coeff.global_integer_offset_y_pix
                ))
                if self.intermediate_results:
                    data_product = np.fft.fft2(data1_rs) * \
                                   np.fft.fft2(data2_rs).conj()
                    cc_image = np.fft.fftshift(np.fft.ifft2(data_product))
                    power = np.log10(cc_image.real)
                    hdu_power = fits.PrimaryHDU(power)
                    hdul_power = fits.HDUList([hdu_power])
                    hdul_power.writeto('power.fits', overwrite=True)
            else:
                rectwv_coeff.global_integer_offset_x_pix = \
                    rinput.global_integer_offset_x_pix
                rectwv_coeff.global_integer_offset_y_pix = \
                    rinput.global_integer_offset_y_pix

            # apply initial rectification and wavelength calibration
            reduced_mos = apply_rectwv_coeff(
                reduced_image,
                rectwv_coeff
            )

            self.logger.info(
                'Refining wavelength calibration (mode={})'.format(
                    rinput.refine_wavecalib_mode
                ))
            # refine RectWaveCoeff object
            rectwv_coeff, expected_catalog_lines = refine_rectwv_coeff(
                reduced_mos,
                rectwv_coeff,
                catlines_all_wave,
                catlines_all_flux,
                rinput.refine_wavecalib_mode,
                list_useful_slitlets,
                save_intermediate_results=self.intermediate_results
            )
            self.save_intermediate_img(expected_catalog_lines,
                                       'expected_catalog_lines.fits')

        # apply rectification and wavelength calibration
        reduced_mos = apply_rectwv_coeff(
            reduced_image,
            rectwv_coeff
        )

        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(rectwv_coeff)
            save_spectral_lines_ds9(rectwv_coeff)

        # compute median spectra employing the useful region of the
        # rectified image
        if self.intermediate_results:
            for imode, outfile in enumerate(['median_spectra_full',
                                             'median_spectra_slitlets',
                                             'median_spectrum_slitlets']):
                median_image = median_slitlets_rectified(
                    reduced_mos, mode=imode
                )
                self.save_intermediate_img(median_image, outfile + '.fits')

        # save results in results directory
        self.logger.info('end rect.+wavecal. reduction of stare spectra')
        result = self.create_result(reduced_mos=reduced_mos,
                                    rectwv_coeff=rectwv_coeff)
        return result

    def set_base_headers(self, hdr):
        newhdr = super(GenerateRectwvCoeff, self).set_base_headers(hdr)
        # Update EXP to 0
        newhdr['EXP'] = 0
        return newhdr


class StareSpectraRectwv(EmirRecipe):
    """Process images in Stare spectra mode applying wavelength calibration

    Note that in this case the wavelength calibration has already been
    determined.

    """
    logger = logging.getLogger(__name__)

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    rectwv_coeff = reqs.RectWaveCoeffRequirement()

    reduced_mos = Result(prods.ProcessedMOS)

    def run(self, rinput):
        self.logger.info('applying existing rect.+wavecal. calibration of '
                         'stare spectra')

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # apply bpm, bias, dark and flat
        reduced_image = basic_processing_with_combination(rinput, flow,
                                                          method=sigmaclip)
        # update header with additional info
        hdr = reduced_image[0].header
        self.set_base_headers(hdr)

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # apply rectification and wavelength calibration
        reduced_mos = apply_rectwv_coeff(
            reduced_image,
            rinput.rectwv_coeff
        )

        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(rinput.rectwv_coeff)
            save_spectral_lines_ds9(rinput.rectwv_coeff)

        # compute median spectra employing the useful region of the
        # rectified image
        if self.intermediate_results:
            for imode, outfile in enumerate(['median_spectra_full',
                                             'median_spectra_slitlets',
                                             'median_spectrum_slitlets']):
                median_image = median_slitlets_rectified(
                    reduced_mos, mode=imode
                )
                self.save_intermediate_img(median_image, outfile + '.fits')

        # save results in results directory
        self.logger.info('end rect.+wavecal. reduction of stare spectra')
        result = self.create_result(reduced_mos=reduced_mos)
        return result

    def set_base_headers(self, hdr):
        newhdr = super(StareSpectraRectwv, self).set_base_headers(hdr)
        # Update EXP to 0
        newhdr['EXP'] = 0
        return newhdr
