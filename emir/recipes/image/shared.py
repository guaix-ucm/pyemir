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

from numina import __version__
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

_logger = logging.getLogger('emir.recipes')

def name_redimensioned_frames(label, step, ext='.fits'):
    dn = '%s_r%s' % (label, ext)
    mn = '%s_mr%s' % (label, ext)
    return dn, mn

def name_object_mask(label, step, ext='.fits'):
    return '%s_mro_i%01d%s' % (label, step, ext)

def name_skybackground(label, step, ext='.fits'):
    dn = '%s_sky_i%01d%s' % (label, step, ext)
    return dn

def name_skybackgroundmask(label, step, ext='.fits'):
    dn = '%s_skymask_i%01d%s' % (label, step, ext)
    return dn

def name_skysub_proc(label, step, ext='.fits'):
    dn = '%s_rfs_i%01d%s' % (label, step, ext)
    return dn

def name_skyflat(label, step, ext='.fits'):
    dn = 'superflat_%s_i%01d%s' % (label, step, ext)
    return dn

def name_skyflat_proc(label, step, ext='.fits'):
    dn = '%s_rf_i%01d%s' % (label, step, ext)
    return dn

def name_segmask(step, ext='.fits'):
    return "check_i%01d%s" % (step, ext)


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

    for idx, frame in enumerate(frames[1:]):
        with pyfits.open(frame) as hdulist:
            wcs = pywcs.WCS(hdulist[0].header)
            pixval = wcs.wcs_sky2pix(skyref, 1)
            result[idx + 1] = pixval[0] - pixref[0]

    return result

class DirectImageCommon(object):
        
    logger = _logger
    
    
    def process(self, obresult, baseshape, amplifiers, subpix=1, 
                store_intermediate=True):
        
        recipe_result = {'products' : []}

        # Convert channels to slices        
        amplifiers = [(slice(*ch0), slice(*ch1)) for ch0, ch1 in amplifiers] 

        if store_intermediate:
            recipe_result['intermediate'] = []
                
        subpixshape = tuple((side * subpix) for side in baseshape)   
        
        # Reference pixel in the center of the frame
        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2)
        
        _logger.info('Computing offsets from WCS information')
        labels = [frame.label for frame in obresult.frames]
        list_of_offsets = offsets_from_wcs(labels, refpix)
        
        # Insert pixel offsets between frames
        for frame, off in zip(obresult.frames, list_of_offsets):
            frame.baselabel = os.path.splitext(frame.label)[0]
            frame.mask = self.parameters['master_bpm']
            # Insert pixel offsets between frames
            frame.pix_offset = off            
            frame.objmask_data = None
            frame.valid_science = True
            _logger.debug('Frame %s, offset=%s', frame.label, off)
        
        # States
        BASIC, PRERED, CHECKRED, FULLRED, COMPLETE = range(5)
        
        state = BASIC
        step = 0
        
        while True:
            if state == BASIC:    
                _logger.info('Basic processing')

                # Basic processing
                
                # FIXME: add this
                bpm = pyfits.getdata(self.parameters['master_bpm'])
                
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

                for frame in obresult.frames:
                    with pyfits.open(frame.label, mode='update') as hdulist:
                            hdulist = basicflow(hdulist)
                  
                          
                state = PRERED
            elif state == PRERED:
                
                _logger.info('Computing relative offsets')
                offsets = [(frame.pix_offset * subpix) for frame in obresult.frames]
                offsets = numpy.round(offsets).astype('int')        
                finalshape, offsetsp = combine_shape(subpixshape, offsets)
                _logger.info('Shape of resized array is %s', finalshape)
                
                # Resizing frames              
                self.resize(obresult.frames, subpixshape, offsetsp, finalshape, 
                            scale=subpix)
                
                # superflat
                _logger.info('Step %d, superflat correction (SF)', step)
                # Compute scale factors (median)           
                self.update_scale_factors(obresult.frames)

                # Create superflat
                superflat = self.compute_superflat(obresult.frames, amplifiers)
            
                # Apply superflat
                self.apply_superflat(obresult.frames, superflat)

                _logger.info('Simple sky correction')
                for frame in obresult.frames:            
                    self.compute_simple_sky(frame)
                
                # Combining the frames
                _logger.info("Step %d, Combining the frames", step)
                
                sf_data = self.combine_frames(obresult.frames)
                      
                _logger.info('Step %d, finished', step)

                state = CHECKRED                
            else:
                break

        hdu = pyfits.PrimaryHDU(sf_data[0])                
        hdr = hdu.header
        hdr.update('NUMXVER', __version__, 'Numina package version')
        hdr.update('NUMRNAM', self.__class__.__name__, 'Numina recipe name')
        hdr.update('NUMRVER', self.__version__, 'Numina recipe version')
        
        hdr.update('FILENAME', 'result.fits')
        hdr.update('IMGTYP', 'TARGET', 'Image type')
        hdr.update('NUMTYP', 'TARGET', 'Data product type')
        
        varhdu = pyfits.ImageHDU(sf_data[1], name='VARIANCE')
        num = pyfits.ImageHDU(sf_data[2], name='MAP')

        result = pyfits.HDUList([hdu, varhdu, num])        
        
        _logger.info("Final frame created")
        recipe_result['products'] = [DataFrame(result), SourcesCatalog()]
        
        return recipe_result    
    
    
    
    
    def compute_simple_sky(self, frame, step=0, save=False):
        
        dst = name_skysub_proc(frame.baselabel, step)
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
            
            _logger.info('Step %d, SC: subtracting sky to frame %s', 
                         step, prev)
            data[frame.valid_region] -= sky

    def combine_frames(self, frames, out=None, step=0):
        _logger.debug('Step %d, opening sky-subtracted frames', step)

        def fits_open(name):
            '''Open FITS with memmap in readonly mode'''
            return pyfits.open(name, mode='readonly', memmap=True)

        frameslll = [fits_open(frame.lastname) for frame in frames if frame.valid_science]
        _logger.debug('Step %d, opening mask frames', step)
        mskslll = [fits_open(frame.resized_mask) for frame in frames if frame.valid_science]
        _logger.debug('Step %d, combining %d frames', step, len(frameslll))
        try:
            extinc = [pow(10, -0.4 * frame.airmass * self.parameters['extinction']) for frame in frames if frame.valid_science]
            data = [i['primary'].data for i in frameslll]
            masks = [i['primary'].data for i in mskslll]
            
            out = quantileclip(data, masks, scales=extinc, dtype='float32', out=out, fclip=0.1)
            
            # saving the three extensions
            pyfits.writeto('result_i%0d.fits' % step, out[0], clobber=True)
            pyfits.writeto('result_var_i%0d.fits' % step, out[1], clobber=True)
            pyfits.writeto('result_npix_i%0d.fits' % step, out[2], clobber=True)
                
            return out
            
        finally:
            _logger.debug('Step %d, closing sky-subtracted frames', step)
            map(lambda x: x.close(), frameslll)
            _logger.debug('Step %d, closing mask frames', step)
            map(lambda x: x.close(), mskslll)
            
    def apply_superflat(self, frames, flatdata, step=0, save=False):
        _logger.info("Step %d, SF: apply superflat", step)

        # Process all frames with the fitted flat
        # FIXME: not sure
        for frame in frames:
            self.correct_superflat(frame, flatdata, step=step, save=save)
        return frames
            
    def correct_superflat(self, frame, fitted, step=0, save=False):
        
        frame.flat_corrected = name_skyflat_proc(frame.baselabel, step)
        
        if save:
            shutil.copyfile(frame.resized_base, frame.flat_corrected)
        else:
            os.rename(frame.resized_base, frame.flat_corrected)
        
        _logger.info("Step %d, SF: apply superflat to image %s", step, frame.flat_corrected)
        with pyfits.open(frame.flat_corrected, mode='update') as hdulist:
            data = hdulist['primary'].data
            datar = data[frame.valid_region]
            data[frame.valid_region] = correct_flatfield(datar, fitted)    
        
        # Copy primary frame extension
        frame.lastname = frame.flat_corrected            
            
    def compute_superflat(self, frames, amplifiers, segmask=None, step=0):
        _logger.info("Step %d, SF: combining the frames without offsets", step)
        try:
            filelist = []
            data = []
            for frame in frames:
                _logger.debug('Step %d, opening resized frame %s', step, frame.resized_base)
                hdulist = pyfits.open(frame.resized_base, memmap=True, mode='readonly')
                filelist.append(hdulist)
                data.append(hdulist['primary'].data[frame.valid_region])

            scales = [frame.median_scale for frame in frames]
            
            masks = None
            if segmask is not None:
                masks = [segmask[frame.valid_region] for frame in frames]
                
            _logger.debug('Step %d, combining %d frames', step, len(data))
            sf_data, _sf_var, sf_num = flatcombine(data, masks, scales=scales, 
                                                    blank=1.0 / scales[0])            
        finally:
            _logger.debug('Step %d, closing resized frames and mask', step)
            for fileh in filelist:               
                fileh.close()            

        # We interpolate holes by channel
        _logger.debug('Step %d, interpolating holes by channel', step)
        for channel in amplifiers:
            mask = (sf_num[channel] == 0)
            if numpy.any(mask):                    
                fixpix2(sf_data[channel], mask, out=sf_data[channel])

        # Normalize, flat has mean = 1
        sf_data /= sf_data.mean()
        
        # Auxiliary data
        sfhdu = pyfits.PrimaryHDU(sf_data)            
        sfhdu.writeto(name_skyflat('comb', step), clobber=True)
        return sf_data
        
    def update_scale_factors(self, frames, step=0):
        _logger.info('Step %d, SF: computing scale factors', step)
        # FIXME: not sure
        for frame in frames:
            region = frame.valid_region
            data = pyfits.getdata(frame.resized_base)[region]
            mask = pyfits.getdata(frame.resized_mask)[region]
            # FIXME: while developing this ::10 is faster, remove later            
            frame.median_scale = numpy.median(data[mask == 0][::10])
            _logger.debug('median value of %s is %f', frame.resized_base, frame.median_scale)
        return frames
        
    def resize(self, frames, baseshape, offsetsp, finalshape, scale=1, step=0):
        _logger.info('Resizing frames and masks')            
        
        for frame, rel_offset in zip(frames, offsetsp):
            region, _ = subarray_match(finalshape, rel_offset, baseshape)
            # Valid region
            frame.valid_region = region
            # Relative offset
            frame.rel_offset = rel_offset
            # names of frame and mask
            framen, maskn = name_redimensioned_frames(frame.baselabel, step)
            frame.resized_base = framen
            frame.resized_mask = maskn
            
            _logger.debug('%s, valid region is %s, relative offset is %s', frame.label, 
                          custom_region_to_str(region), rel_offset)
            self.resize_frame_and_mask(frame, finalshape, framen, maskn, scale)

        return frames

    def resize_frame_and_mask(self, frame, finalshape, framen, maskn, scale):
        _logger.info('Resizing frame %s, subpix x%i', frame.label, scale)
        resize_fits(frame.label, framen, finalshape, frame.valid_region, 
                    scale=scale)

        _logger.info('Resizing mask %s, subpix x%i', frame.label, scale)
        # We don't conserve the sum of the values of the frame here, just
        # expand the mask
        resize_fits(frame.mask, maskn, finalshape, frame.valid_region, 
                    fill=1, scale=scale, conserve=False)
        