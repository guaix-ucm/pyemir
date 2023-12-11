#
# Copyright 2013-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
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

        def flow(x): return x
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
        hdulist = basic_processing_with_combination(
            rinput, flow, method=median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        result = self.create_result(frame=hdulist)

        return result
