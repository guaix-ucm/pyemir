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

Routines shared by image mode recipes

'''
import os
import logging
import shutil

import numpy
import pyfits
import pywcs

from numina.recipes import RecipeBase, Parameter, provides, DataFrame
from numina.flow import SerialFlow, Node
from numina.flow.node import IdNode
from numina.flow.processing import BiasCorrector, FlatFieldCorrector
from numina.flow.processing import DarkCorrector, NonLinearityCorrector, BadPixelCorrector
from numina.array import combine_shape

from numina.image import get_image_shape, resize_fits, custom_region_to_str
from numina.array import combine_shape, correct_flatfield
from numina.array import subarray_match
from numina.array.combine import flatcombine, median, quantileclip

from emir.dataproducts import SourcesCatalog
from emir.dataproducts import create_result

_logger = logging.getLogger('emir.recipes')

def name_redimensioned_images(label, iteration, ext='.fits'):
    dn = '%s_r%s' % (label, ext)
    mn = '%s_mr%s' % (label, ext)
    return dn, mn

def name_object_mask(label, iteration, ext='.fits'):
    return '%s_mro_i%01d%s' % (label, iteration, ext)

def name_skybackground(label, iteration, ext='.fits'):
    dn = '%s_sky_i%01d%s' % (label, iteration, ext)
    return dn

def name_skybackgroundmask(label, iteration, ext='.fits'):
    dn = '%s_skymask_i%01d%s' % (label, iteration, ext)
    return dn

def name_skysub_proc(label, iteration, ext='.fits'):
    dn = '%s_rfs_i%01d%s' % (label, iteration, ext)
    return dn

def name_skyflat(label, iteration, ext='.fits'):
    dn = 'superflat_%s_i%01d%s' % (label, iteration, ext)
    return dn

def name_skyflat_proc(label, iteration, ext='.fits'):
    dn = '%s_rf_i%01d%s' % (label, iteration, ext)
    return dn

def name_segmask(iteration, ext='.fits'):
    return "check_i%01d%s" % (iteration, ext)


def offsets_from_wcs(frames, pixref):
    '''Compute offsets between frames using WCS information.
    
    :parameter frames: sequence of FITS filenames or file descritors
    :parameter pixref: numpy array used as reference pixel
    
    The sky world coordinates are computed on *pixref* using
    the WCS of the first frame in the sequence. Then, the
    pixel coordinates of the reference sky world-coordinates 
    are computed for the rest of the frames.
    
    The results is a numpy array with the difference between the
    computed pixel value and the reference pixel. The first line
    of the array is [0, 0], being the offset from the first image
    to itself. 
    
    '''
    
    result = numpy.zeros((len(frames), pixref.shape[1]))

    with pyfits.open(frames[0]) as hdulist:
        wcs = pywcs.WCS(hdulist[0].header)
        skyref = wcs.wcs_pix2sky(pixref, 1)

    result[0] = pixref[0] - pixref[0]

    for idx, img in enumerate(frames[1:]):
        with pyfits.open(img) as hdulist:
            wcs = pywcs.WCS(hdulist[0].header)
            pixval = wcs.wcs_sky2pix(skyref, 1)
            result[idx + 1] = pixval[0] - pixref[0]

    return result

class RecipeParameters(object):
    pass

class DirectImageCommon(object):
    
    # States
    BASIC, PRERED, CHECKRED, FULLRED, COMPLETE = range(5)
    
    logger = _logger
    
    
    def process(self, obresult, baseshape, amplifiers, subpix=1, 
                store_intermediate=True):
        
        recipe_result = {'products' : []}
        
        if store_intermediate:
            recipe_result['intermediate'] = []
                
        subpixshape = tuple((side * subpix) for side in baseshape)   
        
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
        step = 0
        
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
                offsets = [(frame.pix_offset * subpix) for frame in obresult.frames]
                offsets = numpy.round(offsets).astype('int')        
                finalshape, offsetsp = combine_shape(subpixshape, offsets)
                _logger.info('Shape of resized array is %s', finalshape)
                
                # Resizing images              
                self.resize(obresult.frames, subpixshape, offsetsp, finalshape, 
                            scale=subpix)
                
                # superflat
                _logger.info('Iter %d, superflat correction (SF)', step)
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
                _logger.info("Iter %d, Combining the images", step)
                
                sf_data = self.combine_images(obresult.frames)
                      
                _logger.info('Iter %d, finished', step)

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
    
    
    
    
    def compute_simple_sky(self, frame, itr=0, save=False):
        
        dst = name_skysub_proc(frame.baselabel, itr)
        prev = frame.lastname
        
        if save:
            shutil.copyfile(prev, dst)
        else:
            os.rename(prev, dst)
        
        frame.lastname = dst
        
        with pyfits.open(frame.lastname, mode='update') as hdulist:
            data = hdulist['primary'].data
            valid = data[frame.valid_region]

            if frame.objmask_data is not None:
                _logger.debug('object mask defined')
                msk = frame.objmask_data[valid]
                sky = numpy.median(valid[msk == 0])
            else:
                _logger.debug('object mask empty')
                sky = numpy.median(valid)

            _logger.debug('median sky value is %f', sky)
            frame.median_sky = sky
            
            _logger.info('Iter %d, SC: subtracting sky to frame %s', 
                         itr, prev)
            data[frame.valid_region] -= sky

    def combine_images(self, frames, out=None, itr=0):
        _logger.debug('Iter %d, opening sky-subtracted images', itr)

        def fits_open(name):
            '''Open FITS with memmap in readonly mode'''
            return pyfits.open(name, mode='readonly', memmap=True)

        imgslll = [fits_open(image.lastname) for image in frames if image.valid_science]
        #_logger.debug('Iter %d, opening mask images', itr)
        #mskslll = [fun(image.resized_mask) for image in frames if image.valid_science]
        _logger.debug('Iter %d, combining %d images', itr, len(imgslll))
        try:
            extinc = [pow(10, -0.4 * image.airmass * self.parameters['extinction']) for image in frames if image.valid_science]
            data = [i['primary'].data for i in imgslll]
            #masks = [i['primary'].data for i in mskslll]
            masks = None
            
            out = quantileclip(data, masks, scales=extinc, dtype='float32', out=out, fclip=0.1)
            
            # saving the three extensions
            pyfits.writeto('result_i%0d.fits' % itr, out[0], clobber=True)
            pyfits.writeto('result_var_i%0d.fits' % itr, out[1], clobber=True)
            pyfits.writeto('result_npix_i%0d.fits' % itr, out[2], clobber=True)
                
            return out
            
        finally:
            _logger.debug('Iter %d, closing sky-subtracted images', itr)
            map(lambda x: x.close(), imgslll)
            #_logger.debug('Iter %d, closing mask images', 0)
            #map(lambda x: x.close(), mskslll)
            
    def apply_superflat(self, frames, flatdata, iter=0, save=False):
        _logger.info("Iter %d, SF: apply superflat", iter)

        #copynode = Copy()
        #ffcor = FlatFieldCorrector(flatdata=flatdata)

        # Process all images with the fitted flat
        # FIXME: not sure
        for frame in frames:
            self.correct_superflat(frame, flatdata, iter=iter, save=save)
        return frames
            
    def correct_superflat(self, frame, fitted, iter=0, save=False):
        
        frame.flat_corrected = name_skyflat_proc(frame.baselabel, iter)
        
        if save:
            shutil.copyfile(frame.resized_base, frame.flat_corrected)
        else:
            print frame.resized_base, frame.flat_corrected
            os.rename(frame.resized_base, frame.flat_corrected)
        
        _logger.info("Iter %d, SF: apply superflat to image %s", iter, frame.flat_corrected)
        with pyfits.open(frame.flat_corrected, mode='update') as hdulist:
            data = hdulist['primary'].data
            datar = data[frame.valid_region]
            data[frame.valid_region] = correct_flatfield(datar, fitted)    
        
        # Copy primary image extension
        frame.lastname = frame.flat_corrected            
            
    def compute_superflat(self, frames, amplifiers, segmask=None, iter=0):
        _logger.info("Iter %d, SF: combining the images without offsets", iter)
        try:
            filelist = []
            data = []
            for frame in frames:
                _logger.debug('Iter %d, opening resized frame %s', iter, frame.resized_base)
                hdulist = pyfits.open(frame.resized_base, memmap=True, mode='readonly')
                filelist.append(hdulist)
                data.append(hdulist['primary'].data[frame.valid_region])

            scales = [frame.median_scale for frame in frames]
            
            masks = None
            if segmask is not None:
                masks = [segmask[frame.valid_region] for frame in frames]
                
            _logger.debug('Iter %d, combining %d frames', iter, len(data))
            sf_data, _sf_var, sf_num = flatcombine(data, masks, scales=scales, 
                                                    blank=1.0 / scales[0])            
        finally:
            _logger.debug('Iter %d, closing resized frames and mask', iter)
            for fileh in filelist:               
                fileh.close()            

        # We interpolate holes by channel
        #for channel in amplifiers: 
        #    mask = (sf_num[channel] == 0)
        #    if numpy.any(mask):                    
        #        fixpix2(sf_data[channel], mask, out=sf_data[channel])

        # Normalize, flat has mean = 1
        sf_data /= sf_data.mean()
        
        # Auxilyary data
        sfhdu = pyfits.PrimaryHDU(sf_data)            
        sfhdu.writeto(name_skyflat('comb', iter), clobber=True)
        return sf_data
        
    def update_scale_factors(self, frames, iter=0):
        _logger.info('Iter %d, SF: computing scale factors', iter)
        # FIXME: not sure
        for frame in frames:
            region = frame.valid_region
            data = pyfits.getdata(frame.resized_base)[region]
            #mask = pyfits.getdata(image.resized_mask)[region]
            # FIXME: while developing this ::10 is faster, remove later            
            #image.median_scale = numpy.median(data[mask == 0][::10])
            frame.median_scale = numpy.median(data[::10])
            _logger.debug('median value of %s is %f', frame.resized_base, frame.median_scale)
        return frames
        
    def resize(self, frames, baseshape, offsetsp, finalshape, scale=1, itr=0):
        _logger.info('Resizing images and masks')            
        
        for frame, rel_offset in zip(frames, offsetsp):
            region, _ = subarray_match(finalshape, rel_offset, baseshape)
            # Valid region
            frame.valid_region = region
            # Relative offset
            frame.rel_offset = rel_offset
            # names of frame and mask
            imgn, maskn = name_redimensioned_images(frame.baselabel, itr)
            frame.resized_base = imgn
            frame.resized_mask = maskn
            
            _logger.debug('%s, valid region is %s, relative offset is %s', frame.label, 
                          custom_region_to_str(region), rel_offset)
            self.resize_image_and_mask(frame, finalshape, imgn, maskn, scale)

        return frames

    def resize_image_and_mask(self, image, finalshape, imgn, maskn, scale):
        _logger.info('Resizing image %s, subpix x%i', image.label, scale)
        #resize_fits(image.base, imgn, finalshape, image.region)
        resize_fits(image.label, imgn, finalshape, image.valid_region, scale=scale)

        _logger.info('Resizing mask %s, subpix x%i', image.label, scale)
        #resize_fits(image.label+'mask', maskn, finalshape, image.region, fill=1, scale=scale)