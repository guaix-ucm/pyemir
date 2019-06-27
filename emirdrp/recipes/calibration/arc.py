#
# Copyright 2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Arc calibration Recipe for EMIR
"""

from __future__ import division, print_function

from numina.array.combine import median
from numina.core import Requirement, Result
from numina.types.linescatalog import LinesCatalog

from emirdrp.core.recipe import EmirRecipe
import emirdrp.decorators
from numina.processing.combine import basic_processing_with_combination
from emirdrp.processing.wavecal.apply_rectwv_coeff import apply_rectwv_coeff
from emirdrp.processing.wavecal.rectwv_coeff_from_arc_image \
    import rectwv_coeff_from_arc_image
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import save_four_ds9
import emirdrp.products as prods
import emirdrp.requirements as reqs


class ArcCalibrationRecipe(EmirRecipe):
    """Process arc images applying wavelength calibration"""

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    bound_param = reqs.RefinedBoundaryModelParamRequirement()
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')

    reduced_image = Result(prods.ProcessedImage)
    rectwv_coeff = Result(prods.RectWaveCoeff)
    reduced_55sp = Result(prods.ProcessedMOS)
    reduced_arc = Result(prods.ProcessedMOS)

    @emirdrp.decorators.loginfo
    def run(self, rinput):
        self.logger.info('starting rect.+wavecal. reduction of arc spectra')

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

        # RectWaveCoeff object (with rectification and wavelength calibration
        # coefficients for the particular CSU configuration of the arc image)
        # and HDUList object with the FITS image corresponding to 55 median
        # spectra of each slitlet
        rectwv_coeff, reduced_55sp = rectwv_coeff_from_arc_image(
            reduced_image,
            rinput.bound_param,
            rinput.lines_catalog,
        )

        # generate associated ds9 region files and save them in work directory
        if self.intermediate_results:
            save_four_ds9(rectwv_coeff)

        # apply rectification and wavelength calibration
        reduced_arc = apply_rectwv_coeff(reduced_image, rectwv_coeff)

        # save results in result directory
        self.logger.info('end rect.+wavecal. reduction of arc spectra')
        result = self.create_result(reduced_image=reduced_image,
                                    rectwv_coeff=rectwv_coeff,
                                    reduced_55sp=reduced_55sp,
                                    reduced_arc=reduced_arc)
        return result

    def set_base_headers(self, hdr):
        newhdr = super(ArcCalibrationRecipe, self).set_base_headers(hdr)
        # Update SEC to 0
        newhdr['SEC'] = 0
        return newhdr
