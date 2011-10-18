#
# Copyright 2011 Sergio Pascual
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

# Finite state machine

import logging
import os.path

from numina import RecipeBase, Image, __version__

import naming

__all__ = ['Recipe']

__version__ = '1.0'

_logger = logging.getLogger("emir.recipes")

class FakeRefinePointing(object):
    count = 0

    def __call__(self, images, result):
        _logger.info('Refine pointings')
        _logger.info('Detecting objects from "RESULT"')
        _logger.info('Detecting objects from "IMAGES"')
        if self.count == 0:
            self.count = 1
            return images, True
        else:
            return images, False

class ImageInformation:
    def __init__(self, name, mask):
        self.base = name
        self.label = os.path.splitext(name)[0]
        self.mask = mask
        self.logfile = []

        self.region = None
        self.resized_base = None
        self.resized_mask = None
        self.objmask = None
        self.objmask_data = None       

    def __str__(self):
        return self.base

class Recipe(RecipeBase):
    __requires__ = [Image('master_bpm'),
                    Image('master_bias'),
                    Image('master_dark'),
                    Image('master_flat'),
                    Image('nonlineairy') # FIXME: this is not an image
                    ]
    __provides__ = [Image('science')]

    def __init__(self):
        super(Recipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )
        self.pp = {}
        self.pq = {}

        self.iteration = 0

        self.refine_pointing = FakeRefinePointing()

    def run(self, rb):

        # States
        BASIC, PRERED, CHECKRED, FULLRED, COMPLETE = range(5)

        state = BASIC

        basemask = self.pq['master_bpm']
        images = [ImageInformation(im.path, basemask) for im in rb.images]

        while True:
            if state == BASIC:
                # Basic reduction (dark, nonlin, flat)
                images = self.basic(images)
                state = PRERED
            elif state == PRERED:
                self.iteration += 1        
                images = self.resize(images)
                superflat = self.create_superflat(images)
                images = self.correct_superflat(images, superflat)
                images = self.create_and_correct_sky_1(images)
                result = self.combine(images)

                state = CHECKRED
            elif state == CHECKRED:
                images, recompute = self.refine_pointing(images, result)
                if recompute:
                    _logger.info('Recentering is needed')
                    state = PRERED
                else:
                    _logger.info('Recentering is not needed')
                    images = self.check_photometry(images)
                    state = FULLRED
            elif state == FULLRED:
                self.iteration += 1        
                objmask = self.create_obj_mask(result)
                images = self.update_masks(images, objmask)
    
                superflat = self.create_superflat(images)
                images = self.correct_superflat(images, superflat)
            
                images = self.create_and_correct_sky_2(images)
        
                result = self.combine(images)
        
                if self.iteration >= 4:
                    state = COMPLETE
            else:
                break

        return {'result': {'direct_image': result , 'qa': 1}}

    def create_superflat(self, images):
        name = naming.skyflat('comb', self.iteration)
        _logger.info('Creating superflat %s', name)

    def correct_sky(self):
        _logger.info('Correcting sky %d %d', state, self.iteration)

    def resize(self, images):
        _logger.info('Resizing images to final shape')
        for image in images:
            imgn, maskn = naming.redimensioned_images(image.label, self.iteration)
            image.resized_base = imgn
            image.resized_mask = maskn
            finalshape = None
            self.resize_image_and_mask(image, finalshape, imgn, maskn)
        return images

    def combine(self, images):
        _logger.info('Combining %d images', len(images))
        _logger.info('Result image stored in result_i%0d.fits', self.iteration)
#            pyfits.writeto('result_var_i%0d.fits' % self.iter, out[1], clobber=True)
#            pyfits.writeto('result_npix_i%0d.fits' % self.iter, out[2], clobber=True)

        return None

    def update_masks(self, images, objmask):
        _logger.info('Updating object masks')
        for image in images:
            image.objmask = naming.object_mask(image.label, self.iteration)                            
        return images

    def correct_superflat(self, images, superflat):
        _logger.info('Correcting superflat from images')
        for image in images:
            image.lastname = naming.skyflat_proc(image.label, self.iteration)
            image.flat_corrected = image.lastname
            _logger.debug('Correcting superflat from image %s into %s', image.resized_base, image.flat_corrected)
        return images

    def create_and_correct_sky_1(self, images):
        _logger.info('Removing sky from images')
        for image in images:
            dst = naming.skysub_proc(image.label, self.iteration)
            sky = 0
            image.median_sky = sky
            _logger.debug('Subtracting sky level into image %s', dst)
            image.sky_corrected = dst
            image.lastname = dst
        return images

    def create_and_correct_sky_2(self, images):
        for image in images:
            dst = naming.skysub_proc(image.label, self.iteration)
            _logger.debug('Subtracting sky image into image %s', dst)
            image.sky_corrected = dst
            image.lastname = dst
        return images


    def basicflow(self, image):
        image.logfile.append('Correcting dark')
        _logger.debug('Correcting dark on %s', image)
        image.logfile.append('Correcting non-linearity')
        _logger.debug('Correcting non-linearity on %s', image)
        image.logfile.append('Correcting flat')
        _logger.debug('Correcting flat on %s', image)
        return image

    def basic_processing(self, image, basicflow):
        image = basicflow(image)
        return image

    def basic(self, images):
        _logger.info('Basic processing (dark, non-linearity and flat)')
        for image in images:
            self.basic_processing(image, self.basicflow)
        return images

    def resize_image_and_mask(self, image, finalshape, imgn, maskn):
        _logger.debug('Resizing image %s into %s', image.base, imgn)
        #resize_fits(image.base, imgn, finalshape, image.region)
    
        _logger.debug('Resizing mask %s into %s',image.base, maskn)
        #resize_fits(image.mask, maskn, finalshape, image.region, fill=1)

    def create_obj_mask(self, result):
        _logger.info('Creating object mask')
        return None

    def check_photometry(self, images):
        _logger.info('Checking photometry')
        return images

    def refine_pointing(self, images, result):
        _logger.info('Refine pointings')


        return images, False


