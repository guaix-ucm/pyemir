#
# Copyright 2013-2014 Universidad Complutense de Madrid
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


from numina.core import DataFrame
from numina.core import Requirement, Product, DataProductRequirement
from numina.core.requirements import ObservationResultRequirement

from numina.array.combine import median

from emirdrp.core import EmirRecipe
from emirdrp.products import MasterBias
from emirdrp.products import DataFrameType, MasterIntensityFlat
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement

from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs
from .flows import init_filters_bdf
from .flows import init_filters_bd
from .flows import init_filters_b


_logger = logging.getLogger('numina.recipes.emir')


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
    biasframe = Product(MasterBias)

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
        result = self.create_result(biasframe=DataFrame(hdulist))
        return result


class TestBiasCorrectRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    frame = Product(DataFrameType)

    def run(self, rinput):
        _logger.info('starting simple bias reduction')

        flow = init_filters_b(rinput)
        hdu = basic_processing_with_combination(rinput, flow, method=median)
        hdr = hdu.header
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdulist = fits.HDUList([hdu])

        result = self.create_result(frame=hdulist)
        return result


class TestDarkCorrectRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    frame = Product(DataFrameType)

    def run(self, rinput):
        _logger.info('starting simple dark reduction')

        flow = init_filters_bd(rinput)
        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median)
        hdr = hdulist[0].header
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')

        result = self.create_result(frame=hdulist)

        return result


class TestFlatCorrectRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()

    frame = Product(DataFrameType)

    def run(self, rinput):
        _logger.info('starting simple flat reduction')

        flow = init_filters_bdf(rinput)
        hdulist = basic_processing_with_combination(rinput, flow, method=median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        result = self.create_result(frame=hdulist)

        return result


class TestSkyCorrectRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = DataProductRequirement(MasterIntensityFlat,
                                        'Master Sky calibration'
                                        )

    frame = Product(DataFrameType)

    def run(self, rinput):
        _logger.info('starting simple sky reduction')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow, method=median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        # Update SEC to 0
        hdr['SEC'] = 0
        
        result = self.create_result(frame=hdulist)

        return result
