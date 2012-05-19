#
# Copyright 2008-2012 Universidad Complutense de Madrid
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

'''Recipe for the reduction of imaging mode observations.'''

import logging

from numina.recipes import RecipeBase, Parameter, DataProductParameter
from numina.recipes.requirements import Requirement
from numina.recipes import provides, DataFrame

from emir.dataproducts import MasterBias, MasterDark 
from emir.dataproducts import MasterIntensityFlat, MasterBadPixelMask
from emir.dataproducts import SourcesCatalog, NonLinearityCalibration

from .shared import DirectImageCommon


__author__ = "Sergio Pascual <sergiopr@fis.ucm.es>"

_logger = logging.getLogger("numina.recipes.emir")

@provides(DataFrame, SourcesCatalog)
class DitheredImageRecipe(RecipeBase, DirectImageCommon):
    '''Recipe for the reduction of imaging mode observations.

    Recipe to reduce observations obtained in imaging mode, considering different
    possibilities depending on the size of the offsets between individual images.
    In particular, the following observing modes are considered: stare imaging, nodded
    beamswitched imaging, and dithered imaging. 
    
    A critical piece of information here is a table that clearly specifies which
    images can be labeled as *science*, and which ones as *sky*. Note that some
    images are used both as *science* and *sky* (when the size of the targets are
    small compared to the offsets).
    
    **Observing modes:**
     
     * StareImage
     * Nodded/Beam-switched images
     * Dithered images 
    
    
    **Inputs:**
    
     * Science frames + [Sky Frames]
     * Observing mode name: **stare image**, **nodded beamswitched image**, or **dithered imaging**
     * A table relating each science image with its sky image(s) (TBD if it's in 
       the FITS header and/or in other format)
     * Offsets between them (Offsets must be integer)
     * Master Dark 
     * Bad pixel mask (BPM) 
     * Non-linearity correction polynomials 
     * Master flat (twilight/dome flats)
     * Master background (thermal background, only in K band)
     * Exposure Time (must be the same in all the frames)
     * Airmass for each frame
     * Detector model (gain, RN, lecture mode)
     * Average extinction in the filter
     * Astrometric calibration (TBD)
    
    **Outputs:**
    
     * Image with three extensions: final image scaled to the individual exposure
       time, variance  and exposure time map OR number of images combined (TBD)
    
    **Procedure:**
    
    Images are corrected from dark, non-linearity and flat. Then, an iterative
    process starts:
    
     * Sky is computed from each frame, using the list of sky images of each
       science frame. The objects are avoided using a mask (from the second
       iteration on).
    
     * The relative offsets are the nominal from the telescope. From the second
       iteration on, we refine them using objects of appropriate brightness (not
       too bright, not to faint).
    
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

    __requires__ = [
        DataProductParameter('master_bpm', MasterBadPixelMask, 
                  'Master bad pixel mask'),       
        DataProductParameter('master_bias', MasterBias, 'Master bias image', optional=True),
        DataProductParameter('master_dark', MasterDark, 'Master dark image'),
        DataProductParameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        DataProductParameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        Parameter('extinction', 0.0, 'Mean atmospheric extinction'),
        # FIXME: this parameter is optional 
        Parameter('sources', None, 
                  'List of x, y coordinates to measure FWHM',
                  optional=True),
        Parameter('offsets', None, 'List of pairs of offsets',
                  optional=True),
        Parameter('iterations', 4, 'Iterations of the recipe'),
        Parameter('sky_images', 5, 'Images used to estimate the background before and after current image'),
        Parameter('sky_images_sep_time', 10, 'Maximum separation time between consecutive sky images in minutes'),
        Parameter('check_photometry_levels', [0.5, 0.8], 'Levels to check the flux of the objects'),
        Parameter('check_photometry_actions', ['warn', 'warn', 'default'], 'Actions to take on images'),
        Requirement('instrument.detector.channels', 'List of channels'),
        Requirement('instrument.detector.shape', 'Detector shape'),                    
                    
    ]

    def __init__(self):
        super(DitheredImageRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )
        
    def run(self, obresult):        
        baseshape = self.parameters['instrument.detector.shape']
        amplifiers = self.parameters['instrument.detector.channels']
        offsets = self.parameters['offsets']
        
        
        
        return self.process(obresult, baseshape, amplifiers, 
                            offsets=offsets, subpix=1, 
                            stop_after=DirectImageCommon.FULLRED)
                                                            
