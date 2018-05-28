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


from numina.core import Result
from numina.array.combine import median
import numina.processing as proc

import emirdrp.requirements as reqs
from emirdrp.core.recipe import EmirRecipe
import emirdrp.products as prods
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.processing.wavecal.apply_rectwv_coeff import apply_rectwv_coeff
from emirdrp.processing.wavecal.rectwv_coeff_from_mos_library \
    import rectwv_coeff_from_mos_library
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import save_four_ds9


class StareSpectraRecipe(EmirRecipe):
    """Process images in Stare spectra mode"""

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    master_sky = reqs.SpectralSkyRequirement(optional=True)

    stare = Result(prods.DataFrameType)

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
    master_rectwv = reqs.MasterRectWaveRequirement()
    master_sky = reqs.SpectralSkyRequirement(optional=True)

    reduced_image = Result(prods.DataFrameType)
    stare = Result(prods.DataFrameType)

    def run(self, rinput):
        self.logger.info('starting rect.+wavecal. reduction of stare spectra')

        self.logger.info(rinput.master_rectwv)

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # apply bpm, bias, dark and flat
        reduced_image = basic_processing_with_combination(rinput, flow,
                                                          method=median)
        # update header con additional info
        hdr = reduced_image[0].header
        self.set_base_headers(hdr)

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # RectWaveCoeff object with rectification and wavelength calibration
        # coefficients for the particular CSU configuration
        rectwv_coeff = rectwv_coeff_from_mos_library(
            reduced_image,
            rinput.master_rectwv
        )
        # save as JSON file in work directory
        self.save_structured_as_json(rectwv_coeff, 'rectwv_coeff.json')

        # generate associated ds9 region files and save them in work directory
        if self.intermediate_results:
            save_four_ds9(rectwv_coeff)

        # apply rectification and wavelength calibration
        stare_image = apply_rectwv_coeff(
            reduced_image,
            rectwv_coeff
        )

        if rinput.master_sky:
            # Sky subtraction after rectification
            msky = rinput.master_sky.open()
            sky_corrector = proc.SkyCorrector(
                msky[0].data,
                datamodel=self.datamodel,
                calibid=self.datamodel.get_imgid(msky)
            )

            stare_image = sky_corrector(stare_image)

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
