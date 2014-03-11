#
# Copyright 2011-2014 Universidad Complutense de Madrid
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
from astropy.io import fits
from numina.core import RecipeError
from numina.core import BaseRecipe, RecipeRequirements, DataFrame
from numina.core import Product, DataProductRequirement
from numina.core import define_requirements, define_result
from numina.logger import log_to_history
from numina.array.combine import median
from numina import __version__
from numina.flow import SerialFlow
from numina.flow.node import IdNode
from numina.flow.processing import BiasCorrector
from numina.flow.processing import DarkCorrector, NonLinearityCorrector
# from numina.flow.processing import BadPixelCorrector
from numina.core.requirements import ObservationResultRequirement
from numina.core.requirements import InstrumentConfigurationRequirement

import emir.instrument.channels as allchannels
from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import NonLinearityCalibration, MasterIntensityFlat
from emir.dataproducts import WavelengthCalibration, MasterSpectralFlat
from emir.dataproducts import ChannelLevelStatisticsType
from emir.dataproducts import SlitTransmissionCalibration, ChannelLevelStatistics

__all__ = ['BiasRecipe', 'DarkRecipe', 'IntensityFlatRecipe',
           'SpectralFlatRecipe', 'SlitTransmissionRecipe',
           'WavelengthCalibrationRecipe']

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"

class BiasRecipeRequirements(RecipeRequirements):
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask', optional=True)
    obresult = ObservationResultRequirement()
    insconf = InstrumentConfigurationRequirement()

class BiasRecipeResult(RecipeResult):
    biasframe = Product(MasterBias)
    stats = Product(ChannelLevelStatisticsType)

@define_requirements(BiasRecipeRequirements)
@define_result(BiasRecipeResult)
class BiasRecipe(BaseRecipe):
    '''
    Recipe to process data taken in Bias image Mode.

    Bias images only appear in Simple Readout mode.

    **Observing modes:**
    
     * Bias Image (3.1)   
    
    **Inputs:**
    
    **Outputs:**
    
     * A combined bias frame, with variance extension.
     * Statistics of the final image per channel (mean, median , variance) 
    
    **Procedure:**
    
    The list of images can be readly processed by combining them with a median algorithm.
    
    '''

    def __init__(self):
        super(BiasRecipe, self).__init__(author=_s_author,version="0.1.0")

    # FIXME find a better way of doing this automatically
    # @log_to_history(_logger)
    def run(self, rinput):
        _logger.info('starting bias reduction')

        
        insconf = rinput.insconf.values
        channels_name = insconf['detector']['channels']
        channels = getattr(allchannels, channels_name)
        
        cdata = []
        try:
            for frame in rinput.obresult.frames:
                cdata.append(frame.open())

            _logger.info('stacking %d images using median', len(cdata))
            
            data = median([d['primary'].data for d in cdata], dtype='float32')
            template_header = cdata[0]['PRIMARY'].header
            hdu = fits.PrimaryHDU(data[0], header=template_header)
        finally:
            for hdulist in cdata:
                hdulist.close()

        var2 = numpy.zeros_like(data[0])
        
        def _s_to_f(myslice):
            b = myslice.start
            e = myslice.stop
            return b+1, e
        
        # FIXME: This could be a recarray
        statistics = numpy.empty((len(channels), 7))
        for idx, region in enumerate(channels):
            mean = numpy.mean(data[0][region])
            med = numpy.median(data[0][region])
            var = numpy.var(data[0][region])
            regy, regx = region
            stats = _s_to_f(regy) + _s_to_f(regx) + (mean, med, var)
            statistics[idx, :] = stats 
            var2[region] = var
    
        cls = ChannelLevelStatistics(exposure=0.0, statistics=statistics)
        # update hdu header with
        # reduction keywords
        hdr = hdu.header
        hdr.update('IMGTYP', 'BIAS', 'Image type')
        hdr.update('NUMTYP', 'MASTER_BIAS', 'Data product type')
        hdr.update('NUMXVER', __version__, 'Numina package version')
        hdr.update('NUMRNAM', self.__class__.__name__, 'Numina recipe name')
        hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

        exhdr = fits.Header()
        exhdr.update('extver', 1)
        varhdu = fits.ImageHDU(data[1], name='VARIANCE', header=exhdr)
        exhdr = fits.Header()
        exhdr.update('extver', 2)
        var2hdu = fits.ImageHDU(var2, name='VARIANCE', header=exhdr)
        num = fits.ImageHDU(data[2], name='MAP')

        #hdulist = fits.HDUList([hdu, varhdu, var2hdu, num])
        hdulist = fits.HDUList([hdu])

        _logger.info('bias reduction ended')
            
        result = BiasRecipeResult(biasframe=DataFrame(hdulist),
                                      stats=cls)
        return result
            
class DarkRecipeRequirements(BiasRecipeRequirements):
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)

class DarkRecipeResult(RecipeResult):
    darkframe = Product(MasterDark)
    stats = Product(ChannelLevelStatisticsType)

@define_requirements(DarkRecipeRequirements)
@define_result(DarkRecipeResult)
class DarkRecipe(BaseRecipe):
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

    def __init__(self):
        super(DarkRecipe, self).__init__(author=_s_author, version="0.1.0")
        
    # @log_to_history(_logger)
    def run(self, rinput):
        _logger.info('starting dark reduction')

        output_in_adu_s = True
        
        insconf = rinput.insconf.values
        channels_name = insconf['detector']['channels']
        channels = getattr(allchannels, channels_name)
        
        cdata = []
        expdata = []

        try:
                        
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                exposure = hdulist[0].header['EXPTIME']
                cdata.append(hdulist)
                expdata.append(exposure)

            if not all([exp == expdata[0] for exp in expdata]):
                _logger.error('image with wrong exposure time')
                raise RecipeError('image with wrong exposure time')

            _logger.info('stacking %d images using median', len(cdata))
            
            data = median([d['primary'].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)
            
        finally:
            for hdulist in cdata:
                hdulist.close()

        if rinput.master_bias is not None:
            # load bias
            
            master_bias = rinput.master_bias.open()
            _logger.info('subtracting bias %r', rinput.master_bias)
            # subtrac bias
            data[0] -= master_bias[0].data
            
            #idx = master_bias.index_of(('variance', 1))
            #data[1] += master_bias[idx].data
            
        if output_in_adu_s:
            data[0] /= expdata[0]
            data[1] /= (expdata[0])**2

        var2 = numpy.zeros_like(data[0])

        def _s_to_f(myslice):
            b = myslice.start
            e = myslice.stop
            return b+1, e
        
        statistics = numpy.empty((len(channels), 7))
        for idx, region in enumerate(channels):
            mean = numpy.mean(data[0][region])
            med = numpy.median(data[0][region])
            var = numpy.var(data[0][region])
            regy, regx = region
            stats = _s_to_f(regy) + _s_to_f(regx) + (mean, med, var)
            statistics[idx, :] = stats 
            var2[region] = var
    
        cls = ChannelLevelStatistics(exposure=expdata[0], statistics=statistics)
        # update hdu header with
        # reduction keywords
        hdr = hdu.header
        hdr['NUMXVER'] = ( __version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdr['BUNIT'] = 'ADU/s'
        
        hdr['IMGTYP'] = ('DARK', 'Image type')
        hdr['NUMTYP'] = ('MASTER_DARK', 'Data product type')
        
        exhdr = fits.Header()
        exhdr.update('extver', 1)
        varhdu = fits.ImageHDU(data[1], name='VARIANCE', header=exhdr)
        exhdr = fits.Header()
        exhdr.update('extver', 2)
        var2hdu = fits.ImageHDU(var2, name='VARIANCE', header=exhdr)
        num = fits.ImageHDU(data[2], name='MAP')

        #hdulist = fits.HDUList([hdu, varhdu, var2hdu, num])
        hdulist = fits.HDUList([hdu])

        md = DataFrame(hdulist)

        _logger.info('dark reduction ended')
        result = DarkRecipeResult(darkframe=md,
                                      stats=cls)
        return result


class IntensityFlatRecipeRequirements(DarkRecipeRequirements):
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 'Polynomial for non-linearity correction')

class IntensityFlatRecipeResult(RecipeResult):
    flatframe = Product(MasterIntensityFlat)

@define_requirements(IntensityFlatRecipeRequirements)
@define_result(IntensityFlatRecipeResult)
class IntensityFlatRecipe(BaseRecipe):
    '''Recipe to process data taken in intensity flat-field mode.
        
    Recipe to process intensity flat-fields. The flat-on and flat-off images are
    combined (method?) separately and the subtracted to obtain a thermal subtracted
    flat-field.
    
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
    def __init__(self):
        super(IntensityFlatRecipe, self).__init__(
                        author=_s_author,
                        version="0.1.0"
                )
        
    def run(self, rinput):
        _logger.info('starting flat reduction')
                
        # Basic processing
        if rinput.master_bias:
            _logger.info('loading bias')
            with rinput.master_bias.open() as hdul:
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbiasdata)
        else:
            bias_corrector = IdNode()
            
        with rinput.master_dark.open() as mdark_hdul:
            _logger.info('loading dark')
            mdark = mdark_hdul[0].data
            dark_corrector = DarkCorrector(mdark)

        nl_corrector = NonLinearityCorrector(rinput.nonlinearity)
        
        basicflow = SerialFlow([bias_corrector, dark_corrector, nl_corrector])
       
        cdata = []
        try:
            _logger.info('basic frame reduction')
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                hdulist = basicflow(hdulist)
                cdata.append(hdulist)

            _logger.info('stacking %d images using median', len(cdata))
            
            data = median([d['primary'].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)
            
        finally:
            for hdulist in cdata:
                hdulist.close()

        # update hdu header with
        # reduction keywords
        hdr = hdu.header
        hdr.update('NUMXVER', __version__, 'Numina package version')
        hdr.update('NUMRNAM', self.__class__.__name__, 'Numina recipe name')
        hdr.update('NUMRVER', self.__version__, 'Numina recipe version')
        
        hdr.update('IMGTYP', 'FLAT', 'Image type')
        hdr.update('NUMTYP', 'MASTER_FLAT', 'Data product type')
        
        hdulist = fits.HDUList([hdu])

        md = DataFrame(hdulist)

        result = IntensityFlatRecipeResult(flatframe=md)

        return result
        
        
class SpectralFlatRecipeRequirements(IntensityFlatRecipeRequirements):
    pass

class SpectralFlatRecipeResult(RecipeResult):
    flatframe = Product(MasterSpectralFlat)

@define_requirements(SpectralFlatRecipeRequirements)
@define_result(SpectralFlatRecipeResult)
class SpectralFlatRecipe(BaseRecipe):
    '''Spectral Flatfield Recipe.

    Recipe to process spectral flat-fields. The flat-on and flat-off images are
    combined (method?) separately and the subtracted to obtain a thermal subtracted
    flat-field.

    **Observing modes:**

        * Multislit mask Flat-Field
     
    **Inputs:**

     * A list of lamp-on flats
     * A model of the detector 

    **Outputs:**

     * A combined spectral flat field with with variance extension and quality flag.

    **Procedure:**

     * TBD

    '''

    def __init__(self):
        super(SpectralFlatRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult, rinput):
        return SpectralFlatRecipeResult(flatframe=MasterSpectralFlat(None))

class SlitTransmissionRecipeRequirements(RecipeRequirements):
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask')
    master_bias = DataProductRequirement(MasterBias, 'Master bias image')
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 'Polynomial for non-linearity correction')

class SlitTransmissionRecipeResult(RecipeResult):
    slit = Product(SlitTransmissionCalibration)

@define_requirements(SlitTransmissionRecipeRequirements)
@define_result(SlitTransmissionRecipeResult)
class SlitTransmissionRecipe(BaseRecipe):
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

    def __init__(self):
        super(SlitTransmissionRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    @log_to_history(_logger)
    def run(self, obresult, rinput):
        return SlitTransmissionRecipeResult(slit=SlitTransmissionCalibration())

class WavelengthCalibrationRecipeRequirements(RecipeRequirements):
    master_bpm = DataProductRequirement(MasterBadPixelMask, 'Master bad pixel mask')
    master_bias = DataProductRequirement(MasterBias, 'Master bias image')
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    nonlinearity = DataProductRequirement(NonLinearityCalibration([1.0, 0.0]), 'Polynomial for non-linearity correction')
    master_intensity_ff = DataProductRequirement(MasterIntensityFlat, 'Master intensity flatfield')
    master_spectral_ff = DataProductRequirement(MasterSpectralFlat, 'Master spectral flatfield')
        
class WavelengthCalibrationRecipeResult(RecipeResult):
    cal = Product(WavelengthCalibration)
        
@define_requirements(WavelengthCalibrationRecipeRequirements)
@define_result(WavelengthCalibrationRecipeResult)
class WavelengthCalibrationRecipe(BaseRecipe):
    '''Recipe to calibrate the spectral response.

    **Observing modes:**

        * Wavelength calibration (4.5)
     
    **Inputs:**

     * List of line positions
     * Calibrations upto spectral flatfielding

    **Outputs:**

     * Wavelength calibration structure 

    **Procedure:**

     * TBD
    '''

    def __init__(self):
        super(WavelengthCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    @log_to_history(_logger)
    def run(self, obresult, rinput):
        return WavelengthCalibrationRecipeResult(cal=WavelengthCalibration())


