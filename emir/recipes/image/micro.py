#
# Copyright 2012 Universidad Complutense de Madrid
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
class MicroditheredImageRecipe(RecipeBase):
    '''
    Recipe for the reduction of microdithering imaging.
    
    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing to a number of
    sky positions, with separations of the order of sub arcsecs,
    either by moving the either by nodding the TS, tilting the TS
    M2 or shifting the EMIR DTU, the latter being the most likely
    option. Displacements are of the order of fraction of pixels.
    Images share the large majority of the sky positions so they can
    be coadded. Used for improving the spatial resolution of the
    resulting images and not valid for sky or superflat images.
    

    **Observing modes:**

        * Micro-dithered images
    

    Recipe to reduce observations obtained in imaging mode with microdithering.
    A critical piece of information
    here is a table that clearly specifies which images can be labelled as
    *science*, and which ones as *sky*. Note that some images are used both as
    *science* and *sky* (when the size of the targets are small compared to the
    offsets).

    **Observing modes:**
     * Micro-dithered images 

    **Inputs:**


        * Offsets between them
        * Master Dark 
        * Bad pixel mask (BPM) 
        * Non-linearity correction polynomials 
        * Master flat (twilight/dome flats)
        * Master background (thermal background, only in K band)
        * Detector model (gain, RN, lecture mode)
        * Average extinction in the filter
        * Astrometric calibration (TBD)

    **Outputs:**

     * Image with three extensions: final image scaled to the individual exposure
       time, variance  and exposure time map OR number of images combined (TBD).

    **Procedure:**

    Images are regridded to a integer subdivision of the pixel and then they are
    corrected from dark, non-linearity and flat. It should be desirable that the
    microdithering follows a pattern that can be easily translated to a subdivision
    of the pixel size (by an integer *n* = 2, 3, 4,...) that does not requires a
    too high *n* value. An iterative process starts:

     * Sky is computed from each frame, using the list of sky images of each
       science frame. The objects are avoided using a mask (from the second
       iteration on).

     * The relatiev offsets are the nominal from the telescope. From the second
       iteration on, we refine them using bright objects.

     * We combine the sky-subtracted images, output is: a new image, a variance
       image and a exposure map/number of images used map.

     * An object mask is generated.

     * We recompute the sky map, using the object mask as an additional input. From
       here we iterate (typically 4 times).

     * Finally, the images are corrected from atmospheric extinction and flux
       calibrated.

     * A preliminary astrometric calibration can always be used (using the central
       coordinates of the pointing and the plate scale in the detector). A better
       calibration might be computed using available stars (TBD).

    '''


    logger = _logger

    __requires__ = [
        Parameter('master_bpm', MasterBadPixelMask, 
                  'Master bad pixel mask'),       
        Parameter('master_bias', MasterBias, 'Master bias image', soft=True),
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
        super(MicroditheredImageRecipe, self).__init__(
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

    
