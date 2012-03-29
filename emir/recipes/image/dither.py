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
import os.path
import shutil
import itertools
import math
import operator

import numpy
import pyfits

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection
import numdisplay.zscale
from scipy.spatial import cKDTree as KDTree

import numina.image
from numina.image import get_hdu_shape, resize_fits
from numina.flow import SerialFlow
from numina.image.background import create_background_map
from numina.flow.processing import DarkCorrector, NonLinearityCorrector, BadPixelCorrector
from numina.array import subarray_match
from numina.array import combine_shape, correct_flatfield
from numina.array import fixpix2
from numina.array.combine import flatcombine, median, quantileclip

from numina import __version__
from numina.recipes import Parameter, DataFrame, provides
from numina.recipes import RecipeBase, RecipeError

from numina.util.sextractor import SExtractor
from numina.util.sextractor import open as sopen
import numina.util.sexcatalog as sexcatalog


from emir.dataproducts import create_result, MasterBias, MasterDark 
from emir.dataproducts import MasterIntensityFlat, MasterBadPixelMask
from emir.dataproducts import SourcesCatalog, NonLinearityCalibration
 
from .shared import name_skyflat, name_skyflat_proc
from .shared import name_redimensioned_frames, name_object_mask
from .shared import name_skybackground, name_skybackgroundmask
from .shared import name_skysub_proc, name_segmask
from .shared import DirectImageCommon

#import .instrument.detector as detector

__all__ = ['Recipe']

__author__ = "Sergio Pascual <sergiopr@fis.ucm.es>"

_logger = logging.getLogger("emir.recipes")

#mpl.interactive(True)
mpl.rcParams['toolbar'] = 'None'

def update_sky_related(images, nimages=5):
    
    nox = nimages
    # The first nimages
    for idx, image in enumerate(images[:nox]):
        bw = images[0:2 * nox + 1][:idx]
        fw = images[0:2 * nox + 1][idx + 1:]
        image.sky_related = (bw, fw)

    # Images between nimages and -nimages
    for idx, image in enumerate(images[nox:-nox]):
        bw = images[idx:idx + nox]
        fw = images[idx + nox + 1:idx + 2 * nox + 1]
        image.sky_related = (bw, fw)
 
    # The last nimages
    for idx, image in enumerate(images[-nox:]):
        bw = images[-2 * nox - 1:][0:idx + nox + 1]
        fw = images[-2 * nox - 1:][nox + idx + 2:]
        image.sky_related = (bw, fw)
    
    return images

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
                  soft=True),
        Parameter('offsets', None, 'List of pairs of offsets',
                  soft=True),
        Parameter('iterations', 4, 'Iterations of the recipe'),
        Parameter('sky_images', 5, 'Images used to estimate the background before and after current image'),
        Parameter('sky_images_sep_time', 10, 'Maximum separation time between consecutive sky images in minutes'),
        Parameter('check_photometry_levels', [0.5, 0.8], 'Levels to check the flux of the objects'),
        Parameter('check_photometry_actions', ['warn', 'warn', 'default'], 'Actions to take on images'),                    
                    
    ]

    def __init__(self):
        super(DitheredImageRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )
        
    def run(self, obresult):
                
        baseshape = self.instrument['detectors'][0]
        amplifiers = self.instrument['amplifiers'][0]
        offsets = self.parameters['offsets']
        
        return self.process(obresult, baseshape, amplifiers, 
                            offsets=offsets, subpix=1, stop_after=3)
                                                            
    def _compute_advanced_sky(self, image, objmask):
        '''Create a background map from nearby images.'''
        # Create a copy of the image
        dst = name_skysub_proc(image.baselabel, self.iter)
        shutil.copy(image.lastname, dst)
        image.lastname = dst
        
        # Fraction of julian day
        max_time_sep = self.parameters['sky_images_sep_time'] / 1440.0
        thistime = image.mjd
        
        _logger.info('Iter %d, SC: computing advanced sky for %s', self.iter, image.baselabel)
        desc = []
        data = []
        masks = []
        scales = []

        try:
            idx = 0
            for i in itertools.chain(*image.sky_related):
                time_sep = abs(thistime - i.mjd)
                if time_sep > max_time_sep:
                    _logger.warn('image %s is separated from %s more than %dm', 
                                 i.baselabel, image.baselabel, self.parameters['sky_images_sep_time'])
                    _logger.warn('image %s will not be used', i.baselabel)
                    continue
                filename = i.flat_corrected
                hdulist = pyfits.open(filename, mode='readonly')
                data.append(hdulist['primary'].data[i.valid_region])
                scales.append(numpy.median(data[-1]))
                masks.append(objmask[i.valid_region])
                desc.append(hdulist)
                idx += 1

            _logger.debug('computing background with %d images', len(data))
            sky, _, num = median(data, masks, scales=scales)
            if numpy.any(num == 0):
                # We have pixels without
                # sky background information
                _logger.warn('pixels without sky information in image %s',
                             i.flat_corrected)
                binmask = num == 0
                # FIXME: during development, this is faster
                sky[binmask] = sky[num != 0].mean()
                # To continue we interpolate over the patches
                #fixpix2(sky, binmask, out=sky, iterations=1)
                name = name_skybackgroundmask(image.baselabel, self.iter)
                pyfits.writeto(name, binmask.astype('int16'), clobber=True)
                
            hdulist1 = pyfits.open(image.lastname, mode='update')
            try:
                d = hdulist1['primary'].data[image.valid_region]
                
                # FIXME
                # sky median is 1.0 ?
                sky = sky / numpy.median(sky) * numpy.median(d)
                # FIXME
                self.figure_image(sky, image)                 
                d -= sky
                
                name = name_skybackground(image.baselabel, self.iter)
                pyfits.writeto(name, sky, clobber=True)
                _logger.info('Iter %d, SC: subtracting sky %s to image %s', 
                             self.iter, name, image.lastname)                
            
            finally:
                hdulist1.close()
                                                       
        finally:
            for hdl in desc:
                hdl.close()
                


    

