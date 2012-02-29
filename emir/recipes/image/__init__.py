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

import pyfits
from numina.recipes import RecipeBase, Parameter, provides, DataFrame
from numina.flow import SerialFlow
from numina.flow.node import IdNode
from numina.flow.processing import BiasCorrector
from numina.flow.processing import DarkCorrector, NonLinearityCorrector, BadPixelCorrector

from ...dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from ...dataproducts import MasterIntensityFlat
from ...dataproducts import NonLinearityCalibration
from ...dataproducts import SourcesCatalog
from .base import Recipe as DitheredImageRecipe

__all__ = []

_logger = logging.getLogger('emir.recipes')

@provides(DataFrame, SourcesCatalog)
class StareImageRecipe(RecipeBase):
    '''
    The effect of recording images of the sky in a given pointing
    position of the TS


    **Observing modes:**

        * Stare image
    
    '''

    __requires__ = [
        Parameter('master_bpm', MasterBadPixelMask, 
                  'Master bad pixel mask'),       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        Parameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
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
        
        list_of_wcs = [(img.ra, img.dec) for img in obresult.images]
        list_of_offsets = self.convert_wcs_to_pixels(list_of_wcs)
        
        # Insert pixel offsets between images
        for img, off in zip(obresult.images, list_of_offsets):
            img.pix_offset = off
        
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
                #mflat = pyfits.getdata(self.parameters['master_flat'])
                    
                basicflow = SerialFlow([bias_corrector, 
                                            dark_corrector, 
                                            nl_corrector])

                for img in obresult.images:
                    with pyfits.open(img.label, mode='update') as hdulist:
                            hdulist = basicflow(hdulist)
                            
                state = PRERED
            else:
                break
                
        return {'products': [DataFrame(None), SourcesCatalog()]}
    
    def convert_wcs_to_pixels(self, list_of_wcs):
        '''Convert a list of RA, DEC coordinates to offset in pixels'''
        # FIXME: compute here real values
        return [(0, 0)] * len(list_of_wcs)
    
@provides(DataFrame, SourcesCatalog)
class NBImageRecipe(RecipeBase):
    '''
    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing the TS in cycles
    between two, or more, sky positions. Displacements are larger
    than the EMIR FOV, so the images have no common area. Used
    for sky subtraction


    **Observing modes:**

        * Nodded/Beamswitched images
    
    '''

    __requires__ = [
        Parameter('master_bpm', MasterBadPixelMask, 
                  'Master bad pixel mask'),       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        Parameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        # FIXME: this parameter is optional 
        Parameter('sources', None, 'List of x, y coordinates to measure FWHM')
    ]

    def __init__(self):
        super(NBImageRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DataFrame(None), SourcesCatalog()]}

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

    __requires__ = [
        Parameter('master_bpm', MasterBadPixelMask, 
                  'Master bad pixel mask'),       
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 
                  'Polynomial for non-linearity correction'),
        Parameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        # FIXME: this parameter is optional 
        Parameter('sources', None, 'List of x, y coordinates to measure FWHM')
    ]

    def __init__(self):
        super(MicroditheredImageRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DataFrame(None), SourcesCatalog()]}

@provides(DataFrame, SourcesCatalog)
class MosaicRecipe(RecipeBase):
    '''
    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing to a number of
    sky positions, with separations of the order of the EMIR FOV.
    This command is designed to fully cover a given area on the
    sky, but can also be used to point to a number of sky positions
    on which acquisition is only made at the beginning. Supersky
    frame(s) can be built from the image series.

    **Observing modes:**

        * Mosiac images
    
    '''

    __requires__ = [
        # FIXME: this parameter is optional 
        Parameter('sources', None, 
                  'List of x, y coordinates to measure FWHM')
    ]

    def __init__(self):
        super(MosaicRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DataFrame(None), SourcesCatalog()]}


