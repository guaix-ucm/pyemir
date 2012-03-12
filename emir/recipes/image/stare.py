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

'''
Image mode recipes of EMIR

'''

import logging

#
import os
import shutil
#

import pyfits
import numpy
from numina.recipes import RecipeBase, Parameter, provides, DataFrame
from numina.flow import SerialFlow, Node
from numina.flow.node import IdNode
from numina.flow.processing import BiasCorrector, FlatFieldCorrector
from numina.flow.processing import DarkCorrector, NonLinearityCorrector, BadPixelCorrector
from numina.array import combine_shape

from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import NonLinearityCalibration
from emir.dataproducts import SourcesCatalog
from emir.dataproducts import create_result

from .shared import offsets_from_wcs
from .shared import DirectImageCommon

_logger = logging.getLogger('emir.recipes')

@provides(DataFrame, SourcesCatalog)
class StareImageRecipe(RecipeBase, DirectImageCommon):
    '''
    The effect of recording images of the sky in a given pointing
    position of the TS


    **Observing modes:**

        * Stare image
    
    '''

    logger = _logger

    __requires__ = [
        Parameter('master_bpm', MasterBadPixelMask, 
                  'Master bad pixel mask'),       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        Parameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        Parameter('extinction', 0.0, 'Mean atmospheric extinction'),
        # FIXME: this parameter is optional 
        Parameter('sources', None, 
                  'List of x, y coordinates to measure FWHM',
                  soft=True)
    ]

    def __init__(self):
        super(StareImageRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )
        
    def run(self, obresult):
        
        store_intermediate = True
        
        recipe_result = {'products' : []}
        
        if store_intermediate:
            recipe_result['intermediate'] = []
        
        baseshape = self.instrument['detectors'][0]
        amplifiers = self.instrument['amplifiers'][0]
                
        # Reference pixel in the center of the frame
        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2)
        
        _logger.info('Computing offsets from WCS information')
        labels = [img.label for img in obresult.frames]
        list_of_offsets = offsets_from_wcs(labels, refpix)
        
        # Insert pixel offsets between images
        for img, off in zip(obresult.frames, list_of_offsets):
            img.baselabel = os.path.splitext(img.label)[0]
            # Insert pixel offsets between images
            img.pix_offset = off            
            img.objmask_data = None
            img.valid_science = True
            _logger.debug('Frame %s, offset=%s', img.label, off)
        
        # States
        BASIC, PRERED, CHECKRED, FULLRED, COMPLETE = range(5)
        
        state = BASIC
        
        while True:
            if state == BASIC:    
                _logger.info('Basic processing')

                # Basic processing
                
                # FIXME: add this
                # bpm = pyfits.getdata(self.parameters['master_bpm'])
                if self.parameters['master_bias']:
                    mbias = pyfits.getdata(self.parameters['master_bias'])
                    bias_corrector = BiasCorrector(mbias)
                else:
                    bias_corrector = IdNode()
            
                mdark = pyfits.getdata(self.parameters['master_dark'])
                dark_corrector = DarkCorrector(mdark)
                nl_corrector = NonLinearityCorrector(self.parameters['nonlinearity'])
        
                # FIXME
                mflat = pyfits.getdata(self.parameters['master_intensity_ff'])
                ff_corrector = FlatFieldCorrector(mflat)  
                  
                basicflow = SerialFlow([bias_corrector, 
                                        dark_corrector, 
                                        nl_corrector,
                                        ff_corrector
                                        ])

                for img in obresult.frames:
                    with pyfits.open(img.label, mode='update') as hdulist:
                            hdulist = basicflow(hdulist)
                  
                          
                state = PRERED
            elif state == PRERED:
                
                _logger.info('Computing relative offsets')
                offsets = [frame.pix_offset for frame in obresult.frames]
                offsets = numpy.round(offsets).astype('int')        
                finalshape, offsetsp = combine_shape(baseshape, offsets)
                _logger.info('Shape of resized array is %s', finalshape)
                
                # Resizing images              
                self.resize(obresult.frames, baseshape, offsetsp, finalshape)
                
                # superflat
                _logger.info('Iter %d, superflat correction (SF)', 0)
                # Compute scale factors (median)           
                self.update_scale_factors(obresult.frames)

                # Create superflat
                superflat = self.compute_superflat(obresult.frames, amplifiers)
            
                # Apply superflat
                images_info = self.apply_superflat(obresult.frames, superflat)

                _logger.info('Simple sky correction')
                for image in obresult.frames:            
                    self.compute_simple_sky(image)
                
                # Combining the images
                _logger.info("Iter %d, Combining the images", 0)
                
                sf_data = self.combine_images(obresult.frames)
                      
                _logger.info('Iter %d, finished', 0)

                state = CHECKRED                
            else:
                break

        primary_headers = {'FILENAME': 'result.fits',}
        
        result = create_result(sf_data[0], headers=primary_headers,
                                variance=sf_data[1], 
                                exmap=sf_data[2].astype('int16'))
        
        _logger.info("Final image created")

           
        recipe_result['products'] = [DataFrame(result), SourcesCatalog()]
        
        return recipe_result

    
