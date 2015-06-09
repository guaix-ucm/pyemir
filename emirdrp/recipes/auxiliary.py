#
# Copyright 2011-2015 Universidad Complutense de Madrid
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

'''Auxiliary Recipes for EMIR'''

import logging

import numpy

#from numina.core import BaseRecipeAutoQC
from numina.core import RecipeError
from numina.core import DataFrame#from numina.core import BaseRecipeAutoQC
from numina.core import Product
from numina.logger import log_to_history
from numina.array.combine import median
# from numina.flow.processing import BadPixelCorrector
from numina.core.requirements import ObservationResultRequirement
from numina.core.requirements import InstrumentConfigurationRequirement

import emir.instrument.channels as allchannels
from emir.core import EMIR_BIAS_MODES
from emir.core import gather_info_frames
from emir.core import EmirRecipe
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import WavelengthCalibration, MasterSpectralFlat
from emir.dataproducts import ChannelLevelStatisticsType
from emir.dataproducts import ChannelLevelStatistics
from emir.dataproducts import SlitTransmissionCalibration
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterBadPixelMaskRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
from emir.requirements import MasterSpectralFlatFieldRequirement
from .aiv.flows import init_filters_bdf
from .aiv.flows import init_filters_bd
from .aiv.flows import init_filters_b
from .aiv.flows import basic_processing_with_combination

__all__ = ['BiasRecipe', 'DarkRecipe', 'IntensityFlatRecipe',
           'SpectralFlatRecipe', 'SlitTransmissionRecipe',
           'WavelengthCalibrationRecipe']

_logger = logging.getLogger('numina.recipes.emir')


def _s_to_f(myslice):
    b = myslice.start
    e = myslice.stop
    return b+1, e


class BiasRecipe(EmirRecipe):
    '''
    Recipe to process data taken in Bias image Mode.

    Bias images only appear in Simple Readout mode.

    **Observing modes:**

     * Bias Image (3.1)

    **Inputs:**

    **Outputs:**

     * A combined bias frame, with variance extension.
     * Statistics of the final image per channel (mean, median, variance)

    **Procedure:**

    The list of images can be readly processed by
    combining them with a median algorithm.

    '''
    
    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    insconf = InstrumentConfigurationRequirement()

    biasframe = Product(MasterBias)
    stats = Product(ChannelLevelStatisticsType)

    

    # FIXME find a better way of doing this automatically
    # @log_to_history(_logger)
    def run(self, rinput):
        _logger.info('starting bias reduction')

        insconf = rinput.insconf.values
        channels_name = insconf['detector']['channels']
        channels = getattr(allchannels, channels_name)

        iinfo = gather_info_frames(rinput.obresult.frames)

        if iinfo:
            mode = iinfo[0]['readmode']
            if mode.lower() not in EMIR_BIAS_MODES:
                msg = 'readmode %s, is not a bias mode' % mode
                _logger.error(msg)
                raise RecipeError(msg)

        flow = lambda x: x
        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=False)

        pdata = hdulist[0].data

        statistics = numpy.empty((len(channels), 7))
        for idx, region in enumerate(channels):
            mean = numpy.mean(pdata[region])
            med = numpy.median(pdata[region])
            var = numpy.var(pdata[region])
            regy, regx = region
            stats = _s_to_f(regy) + _s_to_f(regx) + (mean, med, var)
            statistics[idx, :] = stats

        cls = ChannelLevelStatistics(exposure=0.0, statistics=statistics)
        # update hdu header with
        # reduction keywords
        hdr = hdulist[0].header
        hdr['IMGTYP'] = ('BIAS', 'Image type')
        hdr['NUMTYP'] = ('MASTER_BIAS', 'Data product type')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdr['CCDMEAN'] = pdata.mean()

        _logger.info('bias reduction ended')

        result = self.create_result(biasframe=DataFrame(hdulist), stats=cls)
        return result


class DarkRecipe(EmirRecipe):
    '''Recipe to process data taken in Dark current image Mode.

    Recipe to process dark images. The dark images will be combined
    using the median.
    They do have to be of the same exposure time t.

    **Observing mode:**

     * Dark current Image (3.2)

    **Inputs:**

    **Outputs:**

     * A combined dark frame, with variance extension.
    '''

    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    insconf = InstrumentConfigurationRequirement()
    master_bias = MasterBiasRequirement()

    darkframe = Product(MasterDark)
    stats = Product(ChannelLevelStatisticsType)


    # @log_to_history(_logger)
    def run(self, rinput):

        _logger.info('starting dark reduction')
        insconf = rinput.insconf.values
        channels_name = insconf['detector']['channels']
        channels = getattr(allchannels, channels_name)

        flow = init_filters_b(rinput)

        iinfo = gather_info_frames(rinput.obresult.frames)
        ref_exptime = 0.0
        for el in iinfo[1:]:
            if abs(el['texp'] - ref_exptime) > 1e-4:
                _logger.error('image with wrong exposure time')
                raise RecipeError('image with wrong exposure time')

        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=True)

        pdata = hdulist[0].data

        statistics = numpy.empty((len(channels), 7))
        for idx, region in enumerate(channels):
            mean = numpy.mean(pdata[region])
            med = numpy.median(pdata[region])
            var = numpy.var(pdata[region])
            regy, regx = region
            stats = _s_to_f(regy) + _s_to_f(regx) + (mean, med, var)
            statistics[idx, :] = stats
            # var2[region] = var

        cls = ChannelLevelStatistics(exposure=ref_exptime,
                                     statistics=statistics)

        # update hdu header with
        # reduction keywords
        hdr = hdulist[0].header
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')

        hdr['IMGTYP'] = ('DARK', 'Image type')
        hdr['NUMTYP'] = ('MASTER_DARK', 'Data product type')
        hdr['CCDMEAN'] = pdata.mean()

        _logger.info('dark reduction ended')
        result = self.create_result(darkframe=hdulist, stats=cls)
        return result


class IntensityFlatRecipe(EmirRecipe):
    '''Recipe to process data taken in intensity flat-field mode.

    Recipe to process intensity flat-fields. The flat-on and
    flat-off images are combined (method?) separately and the subtracted
    to obtain a thermal subtracted flat-field.

    **Observing modes:**

     * Intensity Flat-Field

    **Inputs:**

      * A master dark frame
      * Non linearity
      * A model of the detector.

    **Outputs:**

     * TBD

    **Procedure:**

     * A combined thermal subtracted flat field, normalized to median 1,
       with with variance extension and quality flag.

    '''

    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    insconf = InstrumentConfigurationRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    flatframe = Product(MasterIntensityFlat)    


    def run(self, rinput):
        _logger.info('starting flat reduction')

        errors = True

        flow = init_filters_bd(rinput)
        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=errors)

        hdr = hdulist[0].header
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')

        hdr['IMGTYP'] = ('FLAT', 'Image type')
        hdr['NUMTYP'] = ('MASTER_FLAT', 'Data product type')

        mm = hdulist[0].data.mean()

        hdulist[0].data /= mm
        if errors:
            hdulist['variance'].data /= (mm*mm)

        result = self.create_result(flatframe=hdulist)

        return result


class SimpleSkyRecipe(EmirRecipe):
    '''Recipe to process data taken in intensity flat-field mode.

    '''

    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    insconf = InstrumentConfigurationRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()

    skyframe = Product(MasterIntensityFlat)


    def run(self, rinput):
        _logger.info('starting sky reduction')

        flow = init_filters_bdf(rinput)

        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=True)

        hdr = hdulist[0].header
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')

        hdr['IMGTYP'] = ('SKY', 'Image type')
        hdr['NUMTYP'] = ('MASTER_SKY', 'Data product type')

        result = self.create_result(skyframe=hdulist)

        return result


class SpectralFlatRecipe(EmirRecipe):
   

    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    insconf = InstrumentConfigurationRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    
    flatframe = Product(MasterSpectralFlat)

    def run(self, obresult, rinput):
        return self.create_result(flatframe=MasterSpectralFlat(None))


class SlitTransmissionRecipe(EmirRecipe):
    '''Recipe to calibrate the slit transmission.

    **Observing modes:**

        * Slit transmission calibration (4.4)

    **Inputs:**

        * A list of uniformly illuminated images of MSM

    **Outputs:**

     * A list of slit transmission functions

    **Procedure:**

     * TBD

    '''

    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    slit = Product(SlitTransmissionCalibration)


    @log_to_history(_logger, 'slit')
    def run(self, obresult, rinput):
        return self.create_result(slit=SlitTransmissionCalibration())



class WavelengthCalibrationRecipe():
    '''Recipe to calibrate the spectral response.

    **Observing modes:**

        * Wavelength calibration (4.5)

    **Inputs:**

     * List of line positions
     * Calibrations up to spectral flatfielding

    **Outputs:**

     * Wavelength calibration structure

    **Procedure:**

     * TBD
    '''

    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_spectral_ff = MasterSpectralFlatFieldRequirement()

    cal = Product(WavelengthCalibration)


    @log_to_history(_logger, 'cal')
    def run(self, obresult, rinput):
        return self.create_result(cal=WavelengthCalibration())
