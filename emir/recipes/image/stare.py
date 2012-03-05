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
import shutil
#

import pyfits
import numpy
from numina.image import get_image_shape, resize_fits
from numina.recipes import RecipeBase, Parameter, provides, DataFrame
from numina.flow import SerialFlow
from numina.flow.node import IdNode
from numina.flow.processing import BiasCorrector, FlatFieldCorrector
from numina.flow.processing import DarkCorrector, NonLinearityCorrector, BadPixelCorrector
from numina.array import combine_shape, correct_flatfield
from numina.array import subarray_match
from numina.array.combine import flatcombine, median, quantileclip

from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import NonLinearityCalibration
from emir.dataproducts import SourcesCatalog
from emir.dataproducts import create_result

from .shared import name_skyflat, name_skyflat_proc
from .shared import name_redimensioned_images, name_skysub_proc
from .shared import offsets_from_wcs

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
        
    def compute_simple_sky(self, image):
        
        dst = name_skysub_proc(image.label, 0)
        prev = image.lastname
        shutil.copy(image.lastname, dst)
        image.lastname = dst
        
        with pyfits.open(image.lastname, mode='update') as hdulist1:
            data = hdulist1['primary'].data
            d = data[image.region]
            
            if image.objmask_data is not None:
                m = image.objmask_data[image.region]
                sky = numpy.median(d[m == 0])
            else:
                _logger.debug('object mask empty')
                sky = numpy.median(d)

            _logger.debug('median sky value is %f', sky)
            image.median_sky = sky
            
            _logger.info('Iter %d, SC: subtracting sky to image %s', 
                         0, prev)
            region = image.region
            data[region] -= sky
                
    def compute_superflat(self, iinfo, segmask, amplifiers):
        _logger.info("Iter %d, SF: combining the images without offsets", 0)
        try:
            filelist = []
            data = []
            for image in iinfo:
                _logger.debug('Iter %d, opening resized image %s', 0, image.resized_base)
                hdulist = pyfits.open(image.resized_base, memmap=True, mode='readonly')
                filelist.append(hdulist)
                data.append(hdulist['primary'].data[image.region])

            scales = [image.median_scale for image in iinfo]

            
            # FIXME: plotting
            #self.figure_median_background(scales)

            masks = None
            if segmask is not None:
                masks = [segmask[image.region] for image in iinfo]
                
            _logger.debug('Iter %d, combining %d images', 0, len(data))
            sf_data, _sf_var, sf_num = flatcombine(data, masks, scales=scales, 
                                                    blank=1.0 / scales[0])            
        finally:
            _logger.debug('Iter %d, closing resized images and mask', 0)
            for fileh in filelist:               
                fileh.close()            

        # We interpolate holes by channel
        #for channel in amplifiers: 
        #    mask = (sf_num[channel] == 0)
        #    if numpy.any(mask):                    
        #        fixpix2(sf_data[channel], mask, out=sf_data[channel])

        # Normalize, flat has mean = 1
        sf_data /= sf_data.mean()
        
        sfhdu = pyfits.PrimaryHDU(sf_data)            
        sfhdu.writeto(name_skyflat('comb', 0), clobber=True)
        return sf_data
        
    def update_scale_factors(self, images_info):

        _logger.info('Iter %d, SF: computing scale factors', 0)
        # FIXME: not sure
        for image in images_info:
            region = image.region
            data = pyfits.getdata(image.resized_base)[region]
            #mask = pyfits.getdata(image.resized_mask)[region]
            # FIXME: while developing this ::10 is faster, remove later            
            #image.median_scale = numpy.median(data[mask == 0][::10])
            image.median_scale = numpy.median(data[::10])
            _logger.debug('median value of %s is %f', image.resized_base, image.median_scale)
        return images_info        
        
    def resize(self, images, baseshape):
        _logger.info('Computing offsets')
        
        offsets = [image.pix_offset for image in images]
        offsets = numpy.round(offsets).astype('int')        
        finalshape, offsetsp = combine_shape(baseshape, offsets)
        _logger.info('Shape of resized array is %s', finalshape)

        _logger.info('Resizing images and masks')            
        
        for image, noffset in zip(images, offsetsp):
            region, _ = subarray_match(finalshape, noffset, baseshape)
            image.region = region
            image.noffset = noffset
            imgn, maskn = name_redimensioned_images(image.label, 0)
            image.resized_base = imgn
            image.resized_mask = maskn
                    
            self.resize_image_and_mask(image, finalshape, imgn, maskn)

        return images

    def resize_image_and_mask(self, image, finalshape, imgn, maskn):
        _logger.info('Resizing image %s', image.label)
        #resize_fits(image.base, imgn, finalshape, image.region)
        resize_fits(image.label, imgn, finalshape, image.region)

        _logger.info('Resizing mask %s', image.label)
        #resize_fits(image.label+'mask', maskn, finalshape, image.region, fill=1)

    def apply_superflat(self, images_info, superflat):
        _logger.info("Iter %d, SF: apply superflat", 0)

        
                    

        # Process all images with the fitted flat
        # FIXME: not sure
        for image in images_info:
            self.correct_superflat(image, superflat)
        return images_info

    def combine_images(self, iinfo, out=None):
        _logger.debug('Iter %d, opening sky-subtracted images', 0)

        def fun(name):
            '''Open FITS with memmap in readonly mode'''
            return pyfits.open(name, mode='readonly', memmap=True)

        imgslll = [fun(image.lastname) for image in iinfo if image.valid_science]
        _logger.debug('Iter %d, opening mask images', 0)
        #mskslll = [fun(image.resized_mask) for image in iinfo if image.valid_science]
        #_logger.debug('Iter %d, combining %d images', 0, len(imgslll))
        try:
            extinc = [pow(10, -0.4 * image.airmass * self.parameters['extinction']) for image in iinfo if image.valid_science]
            data = [i['primary'].data for i in imgslll]
            #masks = [i['primary'].data for i in mskslll]
            masks = None
            if out is not None:
                quantileclip(data, masks, scales=extinc, dtype='float32', out=out, fclip=0.1)
            else:
                out = quantileclip(data, masks, scales=extinc, dtype='float32', fclip=0.1)

            # saving the three extensions
            pyfits.writeto('result_i%0d.fits' % 0, out[0], clobber=True)
            pyfits.writeto('result_var_i%0d.fits' % 0, out[1], clobber=True)
            pyfits.writeto('result_npix_i%0d.fits' % 0, out[2], clobber=True)
                
            return out
            
        finally:
            _logger.debug('Iter %d, closing sky-subtracted images', 0)
            map(lambda x: x.close(), imgslll)
            #_logger.debug('Iter %d, closing mask images', 0)
            #map(lambda x: x.close(), mskslll)


    def correct_superflat(self, image, fitted):
        _logger.info("Iter %d, SF: apply superflat to image %s", 0, image.resized_base)
        with pyfits.open(image.resized_base, mode='readonly') as hdulist:
            data = hdulist['primary'].data[image.region]
            newdata = hdulist['primary'].data.copy()
            newdata[image.region] = correct_flatfield(data, fitted)    
            newheader = hdulist['primary'].header.copy()
        
        # Copy primary image extension
        phdu = pyfits.PrimaryHDU(newdata, newheader)
        image.lastname = name_skyflat_proc(image.label, 0)
        image.flat_corrected = image.lastname
        phdu.writeto(image.lastname, clobber=True)

    def run(self, obresult):
        
        baseshape = self.instrument['detectors'][0]
        amplifiers = self.instrument['amplifiers'][0]
                
        # Reference pixel in the center of the frame
        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2)
        
        labels = [img.label for img in obresult.frames]
        list_of_offsets = offsets_from_wcs(labels, refpix)
        
        # Insert pixel offsets between images
        for img, off in zip(obresult.frames, list_of_offsets):
            img.pix_offset = off
            img.objmask_data = None
            img.valid_science = True
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
                # Resizing images
                print 'resize'                            
                self.resize(obresult.frames, 
                            self.instrument['detectors'][0])
                
                # superflat
                _logger.info('Iter %d, superflat correction (SF)', 0)
                # Compute scale factors (median)           
                self.update_scale_factors(obresult.frames)

                # Create superflat
                superflat = self.compute_superflat(obresult.frames, None, amplifiers)
            
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

           
        return {'products': [DataFrame(result), SourcesCatalog()]}

    
