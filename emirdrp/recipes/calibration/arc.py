#
# Copyright 2018 Universidad Complutense de Madrid
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

"""
Arc calibration Recipe for EMIR
"""

from __future__ import division, print_function

from numina.array.combine import median
from numina.core import Requirement, Parameter, Result
from numina.core.requirements import ObservationResultRequirement
from numina.types.linescatalog import LinesCatalog

from emirdrp.core.recipe import EmirRecipe
import emirdrp.decorators
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.processing.wavecal import apply_rectwv_coeff
from emirdrp.processing.wavecal import rectwv_coeff_from_arc_image
import emirdrp.products as prods
import emirdrp.requirements as reqs


class ArcCalibrationRecipe(EmirRecipe):
    """Process arc images applying wavelength calibration"""

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    bound_param = reqs.RefinedBoundaryModelParamRequirement()
    lines_catalog = Requirement(LinesCatalog, 'Catalog of lines')

    reduced_image = Result(prods.DataFrameType)
    reduced_arc = Result(prods.DataFrameType)
    rectwv_coeff = Result(prods.RectWaveCoeff)

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

        # RectWaveCoeff object with rectification and wavelength calibration
        # coefficients for the particular CSU configuration of the arc image
        rectwv_coeff = rectwv_coeff_from_arc_image(
            reduced_image,
            rinput.bound_param,
            rinput.lines_catalog,
        )

        # apply rectification and wavelength calibration
        reduced_arc = apply_rectwv_coeff(reduced_image, rectwv_coeff)

        # save results in result directory
        self.logger.info('end rect.+wavecal. reduction of arc spectra')
        result = self.create_result(reduced_image=reduced_image,
                                    reduced_arc=reduced_arc,
                                    rectwv_coeff=rectwv_coeff)
        return result

    def set_base_headers(self, hdr):
        newhdr = super(ArcCalibrationRecipe, self).set_base_headers(hdr)
        # Update SEC to 0
        newhdr['SEC'] = 0
        return newhdr
