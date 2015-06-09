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

'''AIV Recipes for EMIR'''

import logging

import numpy
from astropy.io import fits
from scipy.stats import linregress
import matplotlib.pyplot as plt

from numina.core import RecipeError
from numina.core import DataFrame
from numina.core import Requirement, Product, DataProductRequirement
from numina.core.requirements import ObservationResultRequirement

from numina.array.combine import median, mean
from numina import __version__
from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import FlatFieldCorrector, SkyCorrector
from numina.flow.node import IdNode
from numina.flow import SerialFlow

from emir.core import EmirRecipe
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import DataFrameType, MasterIntensityFlat
from emir.dataproducts import DarkCurrentValue, CoordinateList2DType
from emir.core import gather_info
from emir.core import EMIR_BIAS_MODES
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs
from .flows import init_filters_bdf
from .flows import init_filters_bd
from .flows import init_filters_b


_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"


class DarkCurrentRecipe(EmirRecipe):
    '''Recipe to process data taken in Dark Current image Mode.'''

    taskidx = '2.1.03'
    taskname = 'Dark current (+contribution from instrument)'

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()

    darkcurrent = Product(DarkCurrentValue)

    def run(self, rinput):
        _logger.info('starting dark current reduction')

        if rinput.master_bias is not None:
            _logger.debug("master bias=%s", rinput.master_bias)
            master_bias = fits.getdata(rinput.master_bias.filename)
            master_bias_base = master_bias
        else:
            master_bias_base = 0

        values_t = []
        values_d = []
        for frame in rinput.obresult.frames:
            with fits.open(frame.label) as hdulist:
                # FIXME: single images must be corrected to have an uniform
                # exposure time
                texp = hdulist[0].header['exposed']
                corrected = hdulist[0].data - master_bias_base
                corrected_mean = corrected.mean()
                _logger.debug("%s mean=%f exposure=%8.2f", frame.label,
                              corrected_mean, texp)
                values_t.append(texp)
                values_d.append(corrected_mean)

        values_t = numpy.array(values_t)
        values_d = numpy.array(values_d)
        slope, intercept, _r_value, _p_value, _std_err = linregress(values_t,
                                                                    values_d)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Exposure time')
        ax.set_ylabel('Dark current [ADU]')
        ax.plot(values_t, values_d, '-*')
        ax.plot(values_t, slope * values_t + intercept, 'r-')
        fig.savefig('dark-current.png')
        print('slope=', slope, 'intercept=', intercept)

        _logger.info('dark current reduction ended')
        result = self.create_result(darkcurrent=DarkCurrentValue())
        return result


class SimpleBiasRecipe(EmirRecipe):
    '''
    Recipe to process data taken in SimpleBias image Mode.

    Bias images only appear in Simple Readout mode.

    **Outputs:**

     * A combined bias frame, with variance extension.

    **Procedure:**

    The list of images can be readly processed by combining them
    with a median algorithm.
    '''

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