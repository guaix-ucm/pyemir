#
# Copyright 2011-2012 Universidad Complutense de Madrid
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
import time

import numpy
import pyfits
from numina.recipes import RecipeBase, Parameter, provides
from numina.logger import log_to_history

from ..dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from ..dataproducts import NonLinearityCalibration, MasterIntensityFlat
from ..dataproducts import WavelengthCalibration
from ..dataproducts import SlitTransmissionCalibration

__all__ = ['BiasRecipe', 'DarkRecipe', 'FlatRecipe']

_logger = logging.getLogger('emir.recipes')

@provides(MasterBias)
class BiasRecipe(RecipeBase):
    '''
    Recipe to process data taken in Bias image Mode.

    Bias images only appear in Simple Readout mode.

    **Observing modes:**
    
     * Bias Image (3.1)   
    
    **Inputs:**
    
     * A list of bias images
     * A model of the detector (gain, RN)
    
    **Outputs:**
    
     * A combined bias frame, with variance extension and quality flag. 
    
    **Procedure:**
    
    The list of images can be readly processed by combining them with a typical
    sigma-clipping algorithm.
    
    '''
    
    

    __requires__ = [
        Parameter('combine', 'median', 'Combine method'),
        Parameter('resultname', 'result.fits', 'Name of the bias output image'),
    ]

    def __init__(self):
        super(BiasRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    @log_to_history(_logger)
    def run(self, rb):
        _logger.info('starting bias reduction')

        images = rb.images

        cdata = []

        try:
            for image in images:
                hdulist = pyfits.open(image, memmap=True, mode='readonly')
                cdata.append(hdulist)

            _logger.info('stacking images')
            data = numpy.zeros(cdata[0]['PRIMARY'].data.shape, dtype='float32')
            for hdulist in cdata:
                data += hdulist['PRIMARY'].data

            data /= len(cdata)
            data += 2.0

            hdu = pyfits.PrimaryHDU(data, header=cdata[0]['PRIMARY'].header)
    
            # update hdu header with
            # reduction keywords
            hdr = hdu.header
            hdr.update('FILENAME', 'master_bias-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'BIAS', 'Image type')
            hdr.update('NUMTYP', 'MASTER_BIAS', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'BiasRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([hdu])

            _logger.info('bias reduction ended')

            # merge header with HISTORY log
            hdr.ascardlist().extend(history_header.ascardlist())    

            return {'products': [MasterBias(hdulist)]}
        finally:
            for hdulist in cdata:
                hdulist.close()
            
@provides(MasterDark)            
class DarkRecipe(RecipeBase):
    '''Recipe to process data taken in Dark current image Mode.

    Recipe to process dark images. The dark images will be combined 
    weighting with the inverses of the corresponding variance extension. 
    They do not have to be of the same exposure time t, they will be 
    scaled to the same t0 ~ 60s (those with very short exposure time 
    should be avoided).
    
    **Observing mode:**
    
     * Dark current Image (3.2)
    
    **Inputs:**
    
     * A list of dark images 
     * A model of the detector (gain, RN)
    
    **Outputs:**
    
     * A combined dark frame, with variance extension and quality flag. 
    
    **Procedure:**
    
    The process will remove cosmic rays (using a typical sigma-clipping algorithm).
    
    ''' 

    __requires__ = [Parameter('master_bias', MasterBias, 'comment'),
                    Parameter('resultname', 'result.fits', 
                              'Name of the dark output image'),
                    ]

    def __init__(self):
        super(DarkRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )
        
    @log_to_history(_logger)
    def run(self, block):
        _logger.info('starting dark reduction')

        try:
            _logger.info('subtracting bias %s', 
                         str(self.parameters['master_bias']))
            with pyfits.open(self.parameters['master_bias'], mode='readonly') as master_bias:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd['primary'].data
                        data -= master_bias['primary'].data
                

            _logger.info('stacking images from block %d', block.id)

            base = block.images[0]
           
            with pyfits.open(base, memmap=True) as fd:
                data = fd['PRIMARY'].data.copy()
                hdr = fd['PRIMARY'].header
           
            for image in block.images[1:]:
                with pyfits.open(image, memmap=True) as fd:
                    add_data = fd['primary'].data
                    data += add_data

            hdu = pyfits.PrimaryHDU(data, header=hdr)
    
            # update hdu header with
            # reduction keywords
            hdr = hdu.header
            hdr.update('FILENAME', 'master_dark-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'DARK', 'Image type')
            hdr.update('NUMTYP', 'MASTER_DARK', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'DarkRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([hdu])

            _logger.info('dark reduction ended')

            # merge final header with HISTORY log
            hdr.ascardlist().extend(history_header.ascardlist())    

            return {'products': [MasterDark(hdulist)]}
        finally:
            pass

@provides(MasterIntensityFlat)
class IntensityFlatRecipe(RecipeBase):
    '''Recipe to process data taken in intensity flat-field mode.
        
    Recipe to process intensity flat-fields. The flat-on and flat-off images are
    combined (method?) separately and the subtracted to obtain a thermal subtracted
    flat-field.
    
    **Observing modes:**
    
     * Intensity Flat-Field
    
    **Inputs:**
    
      * A list of lamp-on flats
      * A list of lamp-off flats
      * A master dark frame
      * A model of the detector. 
    
    **Outputs:**
    
     * TBD
    
    **Procedure:**
    
     * A combined thermal subtracted flat field, normalized to median 1, 
       with with variance extension and quality flag. 
    
    '''
    __requires__ = [ 
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
    ]
    

    def __init__(self):
        super(IntensityFlatRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )
        
    @log_to_history(_logger)
    def run(self, block):
        _logger.info('starting flat reduction')

        try:
            _logger.info('subtracting bias %s', str(self.parameters['master_bias']))
            with pyfits.open(self.parameters['master_bias'], mode='readonly') as master_bias:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd['primary'].data
                        data -= master_bias['primary'].data
                

            _logger.info('subtracting dark %s', str(self.parameters['master_dark']))
            with pyfits.open(self.parameters['master_dark'], mode='readonly') as master_dark:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd['primary'].data
                        data -= master_dark['primary'].data


            _logger.info('stacking images from block %d', block.id)

            base = block.images[0]
           
            with pyfits.open(base, memmap=True) as fd:
                data = fd['PRIMARY'].data.copy()
                hdr = fd['PRIMARY'].header
           
            for image in block.images[1:]:
                with pyfits.open(image, memmap=True) as fd:
                    add_data = fd['primary'].data
                    data += add_data

            # Normalize flat to mean 1.0
            data[:] = 1.0

            hdu = pyfits.PrimaryHDU(data, header=hdr)

            # update hdu header with
            # reduction keywords
            hdr = hdu.header
            hdr.update('FILENAME', 'master_flat-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'FLAT', 'Image type')
            hdr.update('NUMTYP', 'MASTER_FLAT', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'FlatRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([hdu])

            _logger.info('flat reduction ended')

            # merge final header with HISTORY log
            hdr.ascardlist().extend(history_header.ascardlist())    

            return {'products': [MasterIntensityFlat(hdulist)]}
        finally:
            pass
        
        
@provides(MasterSpectralFlat)
class SpectralFlatRecipe(RecipeBase):
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


    __requires__ = [       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
    ]

    def __init__(self):
        super(Recipe1, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [MasterSpectralFlat(None)]}

@provides(SlitTransmissionCalibration)
class SlitTransmissionRecipe(RecipeBase):
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
    

    __requires__ = [       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
    ]

    def __init__(self):
        super(Recipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    @log_to_history(_logger)
    def run(self, obresult):
        return {'products': [SlitTransmissionCalibration()]}

        
        
@provides(WavelengthCalibration)
class WavelengthCalibrationRecipe(RecipeBase):
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

    __requires__ = [       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        Parameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        Parameter('master_spectral_ff', MasterSpectralFlat, 
                  'Master spectral flatfield'),
    ]

    def __init__(self):
        super(Recipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    @log_to_history(_logger)
    def run(self, obresult):
        return {'products': [WavelengthCalibration()]}
        
