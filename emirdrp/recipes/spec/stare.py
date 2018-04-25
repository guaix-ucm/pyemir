#
# Copyright 2016-2018 Universidad Complutense de Madrid
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
Spectroscopy mode, Stare Spectra
"""


from numina.core import Result
from numina.core.requirements import ObservationResultRequirement
from numina.array.combine import median
import emirdrp.requirements as reqs
from emirdrp.core.recipe import EmirRecipe
import emirdrp.products as prods
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.processing.wavecal import evaluate_rect_wpoly
import emirdrp.decorators

class StareSpectraRecipe(EmirRecipe):
    """Process images in Stare spectra mode"""

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    master_sky = reqs.SpectralSkyRequirement(optional=True)

    stare = Result(prods.DataFrameType)

    @emirdrp.decorators.loginfo
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
    """Process images in Stare spectra mode with wavelength calibration"""

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    master_rectwv = reqs.MasterRectWaveRequirement()
    master_sky = reqs.SpectralSkyRequirement(optional=True)

    reduced_image = Result(prods.DataFrameType)
    stare = Result(prods.DataFrameType)

    @emirdrp.decorators.loginfo
    def run(self, rinput):
        self.logger.info('starting stare spectra reduction')

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
        rectwv_coeff = evaluate_rect_wpoly(
            reduced_image,
            rinput.master_rectwv
        )

        # save results in results directory
        self.logger.info('end stare spectra reduction')
        result = self.create_result(reduced_image=reduced_image,
                                    stare=reduced_image)
        return result

    def set_base_headers(self, hdr):
        newhdr = super(StareSpectraWaveRecipe, self).set_base_headers(hdr)
        # Update SEC to 0
        newhdr['SEC'] = 0
        return newhdr
