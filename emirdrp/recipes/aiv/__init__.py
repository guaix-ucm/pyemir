#
# Copyright 2013-2018 Universidad Complutense de Madrid
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

"""AIV Recipes for EMIR"""

import logging

from astropy.io import fits
from numina.array.combine import median
from numina.core import Result
from numina.core.requirements import ObservationResultRequirement
from numina.processing.combine import basic_processing_with_combination

from emirdrp.core.recipe import EmirRecipe
import emirdrp.products as prods
import emirdrp.requirements as reqs

_logger = logging.getLogger(__name__)


class SimpleBiasRecipe(EmirRecipe):
    """
    Recipe to process data taken in SimpleBias image Mode.

    Bias images only appear in Simple Readout mode.

    **Outputs:**

     * A combined bias frame, with variance extension.

    **Procedure:**

    The list of images can be readly processed by combining them
    with a median algorithm.
    """

    obresult = ObservationResultRequirement()
    biasframe = Result(prods.MasterBias)

    def run(self, rinput):
        _logger.info('starting simple bias reduction')

        flow = lambda x: x
        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=True)

        # update hdu header with
        # reduction keywords
        hdr = hdulist[0].header
        hdr['IMGTYP'] = ('BIAS', 'Image type')
        hdr['NUMTYP'] = ('MASTER_BIAS', 'Data product type')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')

        _logger.info('simple bias reduction ended')

        # qc is QC.UNKNOWN
        result = self.create_result(biasframe=hdulist)
        return result


class TestBiasCorrectRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    frame = Result(prods.ProcessedImage)

    def run(self, rinput):
        _logger.info('starting simple bias reduction')

        flow = self.init_filters(rinput)
        hdu = basic_processing_with_combination(rinput, flow, method=median)
        hdr = hdu[0].header
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdulist = fits.HDUList([hdu])

        result = self.create_result(frame=hdulist)
        return result


class TestRectImageRecipe(EmirRecipe):
    """A Recipe to test GCS handling of rectangular images.

    This appeared as a problem during EMIR first comisioning,
    it was fixed, date 2016-06-21
    """

    obresult = ObservationResultRequirement()
    frame = Result(prods.ProcessedImage)

    def run(self, rinput):
        import numpy
        _logger.info('testing rectangular image')

        data = numpy.zeros((500, 1000), dtype='float32')
        data[200:400, 400:600] = 10000.0
        hdu = fits.PrimaryHDU(data)
        hdulist = fits.HDUList([hdu])
        print("numpy shape of data is", data.shape)

        _logger.info('end testing rectangular image')
        result = self.create_result(frame=hdulist)
        return result


class TestDarkCorrectRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    frame = Result(prods.ProcessedImage)

    def run(self, rinput):
        _logger.info('starting simple dark reduction')

        flow = self.init_filters(rinput)
        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median)
        hdr = hdulist[0].header
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')

        result = self.create_result(frame=hdulist)

        return result


class TestFlatCorrectRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()

    frame = Result(prods.ProcessedImage)

    def run(self, rinput):
        _logger.info('starting simple flat reduction')

        flow = self.init_filters(rinput)
        hdulist = basic_processing_with_combination(rinput, flow, method=median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        result = self.create_result(frame=hdulist)

        return result
