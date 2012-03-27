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
from scipy.spatial import KDTree as KDTree

from numina import __version__
from numina.recipes import RecipeBase, Parameter, provides, DataFrame
from numina.flow import SerialFlow, Node
from numina.flow.node import IdNode
from numina.flow.processing import BiasCorrector, FlatFieldCorrector
from numina.flow.processing import DarkCorrector, NonLinearityCorrector, BadPixelCorrector
from numina.array import combine_shape
from numina.array import fixpix2
from numina.image import resize_fits, custom_region_to_str
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
    
    :parameter frames: sequence of FITS filenames or file descriptors
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

    for idx, frame in enumerate(frames[1:]):
        with pyfits.open(frame) as hdulist:
            wcs = pywcs.WCS(hdulist[0].header)
            pixval = wcs.wcs_sky2pix(skyref, 1)
            result[idx + 1] = pixval[0] - pixref[0]

    return result

def intersection(a, b, scale=1):
    '''Intersection between two segments.'''
    a1, a2 = a
    try:
        b1, b2 = b
    except TypeError:
        b1 = b.start
        b2 = b.stop
        

    if a2 <= b1:
        return None
    if a1 >= b2:
        return None

    # a2 > b1 and a1 < b2

    if a2 <= b2:
        if a1 <= b1:
            return slice(b1 * scale, a2 * scale)
        else:
            return slice(a1 * scale, a2 * scale)
    else:
        if a1 <= b1:
            return slice(b1 * scale, b2 * scale)
        else:
            return slice(a1 * scale, b2 * scale)

def clip_slices(r, region, scale=1):
    '''Intersect slices with a region.''' 
    t = []
    for ch in r:
        a1 = intersection(ch[0], region[0], scale=scale)
        if a1 is None:
            continue
        a2 = intersection(ch[1], region[1], scale=scale)
        if a2 is None:
            continue

        t.append((a1, a2))

    return t

class DirectImageCommon(object):
        
    logger = _logger
    BASIC, PRERED, CHECKRED, FULLRED, COMPLETE = range(5)
    
    
    def process(self, obresult, baseshape, amplifiers, 
                offsets=None, window=None, 
                subpix=1, store_intermediate=True,
                target_is_sky=True, stop_after=PRERED):
        
        # metadata = self.instrument['metadata']
        # FIXME: hardcoded
        metadata = {
         "juliandate": "MJD-OBS", 
         "airmass": "AIRMASS", 
         "detector.mode": "CCDMODE", 
         "filter0": "FILTER", 
         "imagetype": "IMGTYP", 
         "exposure": "EXPTIME"
        }     
        recipe_result = {'products' : []}

        if window is None:
            window = tuple((0, siz) for siz in baseshape)

        if store_intermediate:
            recipe_result['intermediate'] = []
        
        
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
                  
                if stop_after == state:
                    break
                else:
                    state = PRERED
            elif state == PRERED:                
                # Shape of the window
                windowshape = tuple((i[1] - i[0]) for i in window)
                _logger.debug('Shape of window is %s', windowshape)
                # Shape of the scaled window
                subpixshape = tuple((side * subpix) for side in windowshape)
                    
                # Scaled window region
                scalewindow = tuple(slice(*(subpix * i for i in p)) for p in window)
                # Window region
                window = tuple(slice(*p) for p in window)
                
                scaled_amp = clip_slices(amplifiers, window, scale=subpix)

                # Reference pixel in the center of the frame
                refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2).astype('float')
        
                # lists of targets and sky frames
                targetframes = []
                skyframes = []
                
                for frame in obresult.frames:
                    
                    # Getting some metadata from FITS header
                    hdr = pyfits.getheader(frame.label)
                    try:
                        frame.exposure = hdr[str(metadata['exposure'])]
                        #frame.baseshape = get_image_shape(hdr)
                        frame.airmass = hdr[str(metadata['airmass'])]
                        frame.mjd = hdr[str(metadata['juliandate'])]
                    except KeyError as e:
                        raise KeyError("%s in image %s" % (str(e), frame.label))
                    
                    
                    frame.baselabel = os.path.splitext(frame.label)[0]
                    frame.mask = self.parameters['master_bpm']
                    # Insert pixel offsets between frames    
                    frame.objmask_data = None
                    frame.valid_target = False
                    frame.valid_sky = False
                    frame.valid_region = scalewindow
                    if frame.itype == 'TARGET':
                        frame.valid_target = True
                        targetframes.append(frame)
                        if target_is_sky:
                            frame.valid_sky = True
                            skyframes.append(frame)
                    if frame.itype == 'SKY':
                        frame.valid_sky = True
                        skyframes.append(frame)
        
                labels = [frame.label for frame in targetframes]
        
                if offsets is None:
                    _logger.info('Computing offsets from WCS information')
                    
                    list_of_offsets = offsets_from_wcs(labels, refpix)
                else:
                    _logger.info('Using offsets from parameters')
                    list_of_offsets = numpy.asarray(offsets)

                # Insert pixel offsets between frames
                for frame, off in zip(targetframes, list_of_offsets):
                    
                    # Insert pixel offsets between frames
                    frame.pix_offset = off
                    frame.scaled_pix_offset = subpix * off
 
                    _logger.debug('Frame %s, offset=%s, scaled=%s', frame.label, off, subpix * off)

                _logger.info('Computing relative offsets')
                offsets = [(frame.scaled_pix_offset) 
                           for frame in targetframes]
                offsets = numpy.round(offsets).astype('int')        
                finalshape, offsetsp = combine_shape(subpixshape, offsets)
                _logger.info('Shape of resized array is %s', finalshape)
                
                # Resizing target frames        
                self.resize(targetframes, subpixshape, offsetsp, finalshape, 
                            window=window, scale=subpix)
                
                if not target_is_sky:
                    for frame in skyframes:
                        frame.resized_base = frame.label
                        frame.resized_mask = frame.mask    
                
                # superflat
                _logger.info('Step %d, superflat correction (SF)', step)
                # Compute scale factors (median)           
                self.update_scale_factors(obresult.frames)

                # Create superflat
                superflat = self.compute_superflat(obresult.frames, scaled_amp)
            
                # Apply superflat
                self.apply_superflat(obresult.frames, superflat)

                _logger.info('Simple sky correction')
                if target_is_sky:
                    # Each frame is the closest sky frame available
                    
                    for frame in obresult.frames:            
                        self.compute_simple_sky(frame)
                else:
                    self.compute_simple_sky_for_frames(targetframes, skyframes)
                
                # Combining the frames
                _logger.info("Step %d, Combining target frames", step)
                
                sf_data = self.combine_frames(targetframes)
                      
                _logger.info('Step %d, finished', step)

                if stop_after == state:
                    break
                else:
                    state = PRERED
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
    
    def compute_simple_sky_for_frames(self, targetframes, skyframes, 
                            maxsep=5, step=0, save=True):
        
        # build kdtree        
        sarray = numpy.array([frame.mjd for frame in skyframes])
        # shape must be (n, 1)
        sarray = numpy.expand_dims(sarray, axis=1)
        
        # query
        tarray = numpy.array([frame.mjd for frame in targetframes])
        # shape must be (n, 1)
        tarray = numpy.expand_dims(tarray, axis=1)
        
        kdtree = KDTree(sarray)
        
        # 1 / minutes in a day 
        MIN_TO_DAY = 0.000694444
        dis, idxs = kdtree.query(tarray, k=1, 
                                 distance_upper_bound=maxsep * MIN_TO_DAY)
        
        nsky = len(sarray)
        
        for tid, idss in enumerate(idxs):
            try:
                tf = targetframes[tid]
                sf = skyframes[idss]
                self.compute_simple_sky_out(tf, sf, step=step, save=save)
            except IndexError:
                _logger.error('No sky image available for frame %s', tf.lastname)
                raise


    def compute_simple_sky_out(self, frame, skyframe, step=0, save=False):
        _logger.info('Correcting sky in frame %s', frame.lastname)
        _logger.info('with sky computed from frame %s', skyframe.lastname)
        
        if hasattr(skyframe, 'median_sky'):
                sky = skyframe.median_sky
        else:
        
            with pyfits.open(skyframe.lastname, mode='readonly') as hdulist:
                data = hdulist['primary'].data
                valid = data[frame.valid_region]


                if skyframe.objmask_data is not None:
                    _logger.debug('object mask defined')
                    msk = frame.objmask_data[valid]
                    sky = numpy.median(valid[msk == 0])
                else:
                    _logger.debug('object mask empty')
                    sky = numpy.median(valid)

            _logger.debug('median sky value is %f', sky)
            skyframe.median_sky = sky
                
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
            valid -= sky

    def compute_simple_sky(self, frame, step=0, save=False):
        self.compute_simple_sky_out(frame, skyframe=frame, step=0, save=False)

    def combine_frames(self, frames, out=None, step=0):
        _logger.debug('Step %d, opening sky-subtracted frames', step)

        def fits_open(name):
            '''Open FITS with memmap in readonly mode'''
            return pyfits.open(name, mode='readonly', memmap=True)

        frameslll = [fits_open(frame.lastname) for frame in frames if frame.valid_target]
        _logger.debug('Step %d, opening mask frames', step)
        mskslll = [fits_open(frame.resized_mask) for frame in frames if frame.valid_target]
        _logger.debug('Step %d, combining %d frames', step, len(frameslll))
        try:
            extinc = [pow(10, -0.4 * frame.airmass * self.parameters['extinction']) for frame in frames if frame.valid_target]
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
        
    def resize(self, frames, shape, offsetsp, finalshape, window=None, scale=1, step=0):
        _logger.info('Resizing frames and masks')
        for frame, rel_offset in zip(frames, offsetsp):
            if frame.valid_target:
                region, _ = subarray_match(finalshape, rel_offset, shape)
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
                self.resize_frame_and_mask(frame, finalshape, framen, maskn, window, scale)
                
        return frames

    def resize_frame_and_mask(self, frame, finalshape, framen, maskn, window, scale):
        _logger.info('Resizing frame %s, window=%s, subpix=%i', frame.label, 
                     custom_region_to_str(window), scale)
        resize_fits(frame.label, framen, finalshape, frame.valid_region, 
                    window=window, scale=scale)

        _logger.info('Resizing mask %s, subpix x%i', frame.label, scale)
        # We don't conserve the sum of the values of the frame here, just
        # expand the mask
        resize_fits(frame.mask, maskn, finalshape, frame.valid_region, 
                    fill=1, window=window, scale=scale, conserve=False)
        
