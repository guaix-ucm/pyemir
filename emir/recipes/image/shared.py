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

from numina.image import get_image_shape, resize_fits, custom_region_to_str
from numina.array import combine_shape, correct_flatfield
from numina.array import subarray_match
from numina.array.combine import flatcombine, median, quantileclip

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


class DirectImageCommon(object):
    
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