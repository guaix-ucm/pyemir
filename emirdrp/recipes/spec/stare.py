#
# Copyright 2016-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Spectroscopy mode, Stare Spectra
"""


import numpy as np

from numina.core import Result
from numina.core import Requirement, Parameter
from numina.array.combine import median
import numina.processing as proc

import emirdrp.requirements as reqs
from emirdrp.core.recipe import EmirRecipe
import emirdrp.products as prods
from emirdrp.products import RectWaveCoeff
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.processing.wavecal.apply_rectwv_coeff import apply_rectwv_coeff
from emirdrp.processing.wavecal.median_slitlets_rectified \
    import median_slitlets_rectified
from emirdrp.processing.wavecal.rectwv_coeff_from_mos_library \
    import rectwv_coeff_from_mos_library
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import save_four_ds9
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 \
    import save_spectral_lines_ds9
from emirdrp.processing.wavecal.refine_rectwv_coeff import refine_rectwv_coeff

from emirdrp.core import EMIR_MINIMUM_SLITLET_WIDTH_MM
from emirdrp.core import EMIR_MAXIMUM_SLITLET_WIDTH_MM


class StareSpectraRecipe(EmirRecipe):
    """Process images in Stare spectra mode"""

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
                                                    method=median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        # Update SEC to 0
        hdr['SEC'] = 0

        self.logger.info('end stare spectra reduction')
        result = self.create_result(stare=hdulist)
        return result


class StareSpectraWaveRecipe(EmirRecipe):
    """Process images in Stare spectra mode applying wavelength calibration"""

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    master_rectwv = reqs.MasterRectWaveRequirement(optional=True)
    master_sky = reqs.SpectralSkyRequirement(optional=True)
    refine_wavecalib_mode = Parameter(
        0,
        description='Apply wavelength calibration refinement',
        choices=[0, 1, 2, 11, 12]
    )
    minimum_slitlet_width_mm = Parameter(
        float(EMIR_MINIMUM_SLITLET_WIDTH_MM),
        description='Minimum width (mm) for a valid slitlet',
    )
    maximum_slitlet_width_mm = Parameter(
        float(EMIR_MAXIMUM_SLITLET_WIDTH_MM),
        description='Maximum width (mm) for a valid slitlet',
    )
    global_integer_offset_x_pix = Parameter(
        0,
        description='Global offset (pixels) in wavelength direction (integer)',
    )
    global_integer_offset_y_pix = Parameter(
        0,
        description='Global offset (pixels) in spatial direction (integer)',
    )

    reduced_image = Result(prods.ProcessedImage)
    stare = Result(prods.ProcessedMOS)

    def run(self, rinput):
        self.logger.info('starting rect.+wavecal. reduction of stare spectra')

        self.logger.info(rinput.master_rectwv)
        self.logger.info('Wavelength calibration refinement mode: {}'.format(
            rinput.refine_wavecalib_mode))
        self.logger.info('Minimum slitlet width (mm)............: {}'.format(
            rinput.minimum_slitlet_width_mm))
        self.logger.info('Maximum slitlet width (mm)............: {}'.format(
            rinput.maximum_slitlet_width_mm))
        self.logger.info('Global offset X direction (pixels)....: {}'.format(
            rinput.global_integer_offset_x_pix))
        self.logger.info('Global offset Y direction (pixels)....: {}'.format(
            rinput.global_integer_offset_y_pix))

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # apply bpm, bias, dark and flat
        reduced_image = basic_processing_with_combination(rinput, flow,
                                                          method=median)
        # update header with additional info
        hdr = reduced_image[0].header
        self.set_base_headers(hdr)

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        stare_image = reduced_image
        if rinput.master_rectwv:
            # RectWaveCoeff object with rectification and wavelength
            # calibration coefficients for the particular CSU configuration
            rectwv_coeff = rectwv_coeff_from_mos_library(
                reduced_image,
                rinput.master_rectwv
            )

            # set global offsets
            rectwv_coeff.global_integer_offset_x_pix = \
                rinput.global_integer_offset_x_pix
            rectwv_coeff.global_integer_offset_y_pix = \
                rinput.global_integer_offset_y_pix

            # apply rectification and wavelength calibration
            stare_image = apply_rectwv_coeff(
                reduced_image,
                rectwv_coeff
            )

            # wavelength calibration refinement
            # 0 -> no refinement
            # 1 -> apply global offset to all the slitlets (using ARC lines)
            # 2 -> apply individual offset to each slitlet (using ARC lines)
            # 11 -> apply global offset to all the slitlets (using OH lines)
            # 12 -> apply individual offset to each slitlet (using OH lines)
            if rinput.refine_wavecalib_mode != 0:
                self.logger.info(
                    'Refining wavelength calibration (mode={})'.format(
                        rinput.refine_wavecalib_mode
                    ))
                # refine RectWaveCoeff object
                rectwv_coeff, expected_catalog_lines = refine_rectwv_coeff(
                    stare_image,
                    rectwv_coeff,
                    rinput.refine_wavecalib_mode,
                    rinput.minimum_slitlet_width_mm,
                    rinput.maximum_slitlet_width_mm,
                    save_intermediate_results=self.intermediate_results
                )
                self.save_intermediate_img(expected_catalog_lines,
                                           'expected_catalog_lines.fits')
                # re-apply rectification and wavelength calibration
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
                self.logger.warning("sky and current image don't have the same shape")
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
        self.logger.info('end rect.+wavecal. reduction of stare spectra')
        result = self.create_result(reduced_image=reduced_image,
                                    stare=stare_image)
        return result

    def set_base_headers(self, hdr):
        newhdr = super(StareSpectraWaveRecipe, self).set_base_headers(hdr)
        # Update SEC to 0
        newhdr['SEC'] = 0
        return newhdr


class StareSpectraApplyWaveRecipe(EmirRecipe):
    """Process images in Stare spectra mode applying wavelength calibration

    Note that in this case the wavelength calibration has already been
    determined.

    """

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    rectwv_coeff_json = Parameter(
        'undefined_rectwv_coeff.json',
        description='JSON file with rect.+wavecal. coefficients',
        optional=False
    )

    reduced_image = Result(prods.ProcessedImage)
    stare = Result(prods.ProcessedMOS)

    def run(self, rinput):
        self.logger.info('applying existing rect.+wavecal. calibration of '
                         'stare spectra')

        self.logger.info('JSON file with rect.+wavecal. coeff.: {}'.format(
            rinput.rectwv_coeff_json))
        # generate RectWaveCoeff object
        rectwv_coeff = RectWaveCoeff._datatype_load(
            '../data/'+rinput.rectwv_coeff_json)

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # apply bpm, bias, dark and flat
        reduced_image = basic_processing_with_combination(rinput, flow,
                                                          method=median)
        # update header with additional info
        hdr = reduced_image[0].header
        self.set_base_headers(hdr)

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # apply rectification and wavelength calibration
        stare_image = apply_rectwv_coeff(
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
                    stare_image, mode=imode
                )
                self.save_intermediate_img(median_image, outfile + '.fits')

        # save results in results directory
        self.logger.info('end rect.+wavecal. reduction of stare spectra')
        result = self.create_result(reduced_image=reduced_image,
                                    stare=stare_image)
        return result

    def set_base_headers(self, hdr):
        newhdr = super(StareSpectraApplyWaveRecipe, self).set_base_headers(hdr)
        # Update SEC to 0
        newhdr['SEC'] = 0
        return newhdr
