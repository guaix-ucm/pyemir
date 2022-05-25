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
from numina.processing.combine import basic_processing_with_combination

from emirdrp.core.recipe import EmirRecipe
import emirdrp.requirements as reqs
import emirdrp.products as prods
from emirdrp.processing.wavecal.apply_rectwv_coeff import apply_rectwv_coeff
from emirdrp.processing.wavecal.rectwv_coeff_from_mos_library \
    import rectwv_coeff_from_mos_library
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import save_four_ds9


class SkySpecRecipe(EmirRecipe):
    """Recipe to process data taken in spectral sky mode.

    """

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    master_rectwv = reqs.MasterRectWaveRequirement()
    skyspec = Result(prods.SkySpectrum)
    reduced_image = Result(prods.ProcessedMOS)

    def run(self, rinput):
        self.logger.info('starting spectral sky reduction')

        flow = self.init_filters(rinput)

        reduced_image = basic_processing_with_combination(
            rinput, flow,
            method=median,
            errors=True
        )

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
        skyspec = apply_rectwv_coeff(
            reduced_image,
            rectwv_coeff
        )

        self.logger.info('end sky spectral reduction')

        result = self.create_result(
            reduced_image=reduced_image,
            skyspec=skyspec
        )

        return result
