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

VALID_SCIENCE = 1
VALID_SKY = 2

class ImageInformation(object):
    '''Selected metadata from an image.
    
    Selected metadata from an image. It stores also processing status'''
    def __init__(self):
        self.baselabel = ''
        self.airmass = 0
        self.mjd = 0
        self.exposure = 0 # Exposure time
        self.flags = VALID_SCIENCE | VALID_SKY
        self.valid_science = True
        self.valid_sky = True

@provides(DataFrame, SourcesCatalog)
class DitheredImageRecipe(RecipeBase, DirectImageCommon):
    '''
    The effect of recording images of the sky in a given pointing
    position of the TS


    **Observing modes:**

        * Stare image
    
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
        
@provides(DataFrame)
@provides(SourcesCatalog)
class DitheredImageRecipe2(RecipeBase):
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
        Parameter('extinction', 0.0, 'Mean atmospheric extinction'),
        Parameter('master_bias', MasterBias, 'Master bias image'),
        Parameter('master_dark', MasterDark, 'Master dark image'),
        Parameter('master_bpm', MasterBadPixelMask, 'Master bad pixel mask'),
        Parameter('master_intensity_ff', MasterIntensityFlat, 
                  'Master intensity flatfield'),
        Parameter('nonlinearity', NonLinearityCalibration([1.0, 0.0]), 'Polynomial for non-linearity correction'),
        Parameter('iterations', 4, 'Iterations of the recipe'),
        Parameter('sky_images', 5, 'Images used to estimate the background before and after current image'),
        Parameter('sky_images_sep_time', 10, 'Maximum separation time between consecutive sky images in minutes'),
        Parameter('resultname', 'result.fits', 'Name of the output image'),
        Parameter('check_photometry_levels', [0.5, 0.8], 'Levels to check the flux of the objects'),
        Parameter('check_photometry_actions', ['warn', 'warn', 'default'], 'Actions to take on images'),
    ]    
    
    def __init__(self):
        super(DitheredImageRecipe, self).__init__(author=__author__, version="0.1.0")
        
        self.iter = 0
        self._figure = plt.figure(facecolor='white')
        self._figure.canvas.set_window_title('Recipe Plots')
        self._figure.canvas.draw()
                            
    def basic_processing(self, image, flow):
        with pyfits.open(image.base) as hdulist:
            hdu = hdulist['primary']         

            # Processing
            hdu = flow(hdu)
            
            hdulist.writeto(os.path.basename(image.base), clobber=True)
            
    def compute_advanced_sky2(self, image, objmask):
        '''Create a background map of the image and subtract it.'''
        # Create a copy of the image
        dst = name_skysub_proc(image.baselabel, self.iter)
        shutil.copy(image.lastname, dst)
        image.lastname = dst
        
        data = pyfits.getdata(image.lastname)
        sky, _ = create_background_map(data[image.valid_region], 128, 128)
        
        _logger.info('Computing advanced sky for %s', image.baselabel)            
        
        with pyfits.open(image.lastname, mode='update') as hdulist:
            d = hdulist['primary'].data[image.valid_region]
        
            _logger.info('Iter %d, SC: subtracting sky to image %s', 
                     self.iter, image.baselabel)
            name = name_skybackground(image.baselabel, self.iter)
            pyfits.writeto(name, sky, clobber=True)
            d -= sky            
            
    def compute_advanced_sky(self, image, objmask):
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
       
    def figures_after_combine(self, sf_data):
       
         # FIXME, more plots
        def truncated(array, frac=0.1):
            '''Dirty truncated mean'''
            nf = int(array.size * frac)
            array.sort()
            new = array[nf:-nf]
            return new.mean(), new.std()
            
        ndata = sf_data[2].astype('int')                        
        data = sf_data[0]
            
        nimages = ndata.max()

        rnimage = range(1, nimages + 1)
        rmean = rnimage[:]
        rstd = rnimage[:]
            
        for pix in rnimage:
            rmean[pix - 1], rstd[pix - 1] = truncated(data[ndata == pix])                
            
        avg_rms = self.figure_check_combination(rnimage, rmean, rstd)
                        
        # Fake sky error image
        self.figure_simple_image(sf_data[2], title='Number of images combined')

        # Create fake error image
        fake = numpy.where(sf_data[2] > 0, numpy.random.normal(avg_rms / numpy.sqrt(sf_data[2])), 0.0)
        self.figure_simple_image(fake, title='Fake sky error image')
        # store fake image
        pyfits.writeto('fake_sky_rms_i%0d.fits' % step, fake, clobber=True)



         
    def figure_check_combination(self, rnimage, rmean, rstd):            
        self._figure.clf()
        self._figure.subplots_adjust(hspace=0.001)
        
        ax1 = self._figure.add_subplot(3,1,1)
        pred = [rstd[-1] * math.sqrt(rnimage[-1] / float(npix)) for npix in rnimage]
        ax1.plot(rnimage, rstd, 'g*', rnimage, pred, 'y-')
        ax1.set_title("")
        ax1.set_ylabel('measured sky rms')
        
        ax2 = self._figure.add_subplot(3,1,2, sharex=ax1)
        pred = [val * math.sqrt(npix) for val, npix in zip(rstd, rnimage)]
        avg_rms = sum(pred) / len(pred)
        ax2.plot(rnimage, pred, 'r*', [rnimage[0], rnimage[-1]], [avg_rms,avg_rms])
        ax2.set_ylabel('scaled sky rms')

        ax3 = self._figure.add_subplot(3,1,3, sharex=ax1)
        ax3.plot(rnimage, rmean, 'b*')
        ax3.set_ylabel('mean sky')
        ax3.set_xlabel('number of frames per pixel')

        xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
        mpl.artist.setp(xticklabels, visible=False)
        self._figure.canvas.draw()
        self._figure.savefig('figure-check-combination_i%01d.png' % step)
        return avg_rms         

    def figure_simple_image(self, data, title=None):
        self._figure.clf()
        ax = self._figure.add_subplot(111)
        cmap = mpl.cm.get_cmap('gray')
        norm = mpl.colors.LogNorm()
        if title is not None:
            ax.set_title(title)
                          
        ax.set_xlabel('X')
        ax.set_ylabel('Y')            
        ax.imshow(data, cmap=cmap, norm=norm)                                
        self._figure.canvas.draw()

    def figure_final_before_s(self, data):
        self._figure.clf()
        ax = self._figure.add_subplot(111)
        cmap = mpl.cm.get_cmap('gray')
        norm = mpl.colors.LogNorm()
        ax.set_title('Result image')              
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        z1, z2 = numdisplay.zscale.zscale(data)
        ax.imshow(data, cmap=cmap, clim=(z1, z2))                                
        self._figure.canvas.draw()

    def figure_fwhm_histogram(self, fwhms):
        self._figure.clf()
        ax = self._figure.add_subplot(111)
        ax.set_title('FWHM of objects')
        ax.hist(fwhms, 50, normed=1, facecolor='g', alpha=0.75)
        self._figure.canvas.draw()
        self._figure.savefig('figure-fwhm-histogram_i%01d.png' % step)
                   
    def figure_init(self, shape):
        self._figure.clf()
        ax = self._figure.add_subplot(111)
        cmap = mpl.cm.get_cmap('gray')
        norm = mpl.colors.LogNorm()
        ax.imshow(numpy.zeros(shape), cmap=cmap, norm=norm)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        self._figure.canvas.draw()      
                 
    def figure_image(self, thedata, image):        
        ax = self._figure.gca()
        image_axes, = ax.get_images()
        image_axes.set_data(thedata)
        z1, z2 = numdisplay.zscale.zscale(thedata)
        image_axes.set_clim(z1, z2)
        clim = image_axes.get_clim()
        ax.set_title('%s, bg=%g fg=%g, linscale' % (image.lastname, clim[0], clim[1]))        
        self._figure.canvas.draw()
        
    def figure_median_background(self, scales):
        # FIXME: plotting
        self._figure.clf()
        ax = self._figure.add_subplot(1,1,1) 
        ax.plot(scales, 'r*')
        ax.set_title("")
        ax.set_xlabel('Image number')
        ax.set_ylabel('Median')
        self._figure.canvas.draw()
        self._figure.savefig('figure-median-sky-background_i%01d.png' % step)

    def compute_superflat(self, iinfo, segmask, amplifiers):
        _logger.info("Iter %d, SF: combining the images without offsets", step)
        try:
            filelist = []
            data = []
            for image in iinfo:
                _logger.debug('Iter %d, opening resized image %s', step, image.resized_base)
                hdulist = pyfits.open(image.resized_base, memmap=True, mode='readonly')
                filelist.append(hdulist)
                data.append(hdulist['primary'].data[image.valid_region])

            scales = [image.median_scale for image in iinfo]

            
            # FIXME: plotting
            self.figure_median_background(scales)

            masks = None
            if segmask is not None:
                masks = [segmask[image.valid_region] for image in iinfo]
                
            _logger.debug('Iter %d, combining %d images', step, len(data))
            sf_data, _sf_var, sf_num = flatcombine(data, masks, scales=scales, 
                                                    blank=1.0 / scales[0])            
        finally:
            _logger.debug('Iter %d, closing resized images and mask', step)
            for fileh in filelist:               
                fileh.close()            

        # We interpolate holes by channel
        for channel in amplifiers: 
            mask = (sf_num[channel] == 0)
            if numpy.any(mask):                    
                fixpix2(sf_data[channel], mask, out=sf_data[channel])

        # Normalize, flat has mean = 1
        sf_data /= sf_data.mean()
        
        sfhdu = pyfits.PrimaryHDU(sf_data)            
        sfhdu.writeto(name_skyflat('comb', step), clobber=True)
        return sf_data

    def correct_superflat(self, image, fitted):
        _logger.info("Iter %d, SF: apply superflat to image %s", step, image.resized_base)
        hdulist = pyfits.open(image.resized_base, mode='readonly')
        data = hdulist['primary'].data[image.valid_region]
        newdata = hdulist['primary'].data.copy()
        newdata[image.valid_region] = correct_flatfield(data, fitted)

                
        newheader = hdulist['primary'].header.copy()
        hdulist.close()
        phdu = pyfits.PrimaryHDU(newdata, newheader)
        image.lastname = name_skyflat_proc(image.baselabel, step)
        image.flat_corrected = image.lastname
        phdu.writeto(image.lastname, clobber=True)
        
        # FIXME: plotting
        try:
            self.figure_image(newdata[image.valid_region], image)
        except ValueError:
            _logger.warning('Problem plotting %s', image.lastname)
        
    def check_photometry_plot(self, vals, errors, levels, nsigma):
        x = range(len(errors))
        self._figure.clf()
        ax = self._figure.add_subplot(111)
        ax.set_title('Relative flux of brightest object')
        for v,c in zip(vals, ['b', 'r', 'g', 'y']):
            ax.scatter(v[0], v[1], c=c)
            w = errors[v[0]]
            ax.errorbar(v[0], v[1], yerr=w, fmt=None, c=c)
            

        ax.plot([x[0], x[-1]], [1, 1], 'r--')
        ax.plot([x[0], x[-1]], [1 - nsigma, 1 - nsigma], 'b--')
        for f in levels:
            ax.plot([x[0], x[-1]], [f, f], 'g--')
            
        self._figure.canvas.draw()
        self._figure.savefig('figure-relative-flux_i%01d.png' % step)
        
    def check_photometry(self, images_info, sf_data, seeing_fwhm):
        # Check photometry of few objects
        weigthmap = 'weights4rms.fits'
        
        wmap = numpy.zeros_like(sf_data[0])
        
        # Center of the image
        border = 300
        wmap[border:-border, border:-border] = 1                    
        pyfits.writeto(weigthmap, wmap, clobber=True)
        
        basename = 'result_i%0d.fits' % (step)
        sex = SExtractor()
        sex.config['VERBOSE_TYPE'] = 'QUIET'
        sex.config['PIXEL_SCALE'] = 1
        sex.config['BACK_TYPE'] = 'AUTO' 
        if seeing_fwhm is not None:
            sex.config['SEEING_FWHM'] = seeing_fwhm * sex.config['PIXEL_SCALE']
        sex.config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
        sex.config['WEIGHT_IMAGE'] = weigthmap
        
        sex.config['PARAMETERS_LIST'].append('FLUX_BEST')
        sex.config['PARAMETERS_LIST'].append('FLUXERR_BEST')
        sex.config['PARAMETERS_LIST'].append('FWHM_IMAGE')
        sex.config['PARAMETERS_LIST'].append('CLASS_STAR')
        
        sex.config['CATALOG_NAME'] = 'master-catalogue-i%01d.cat' % step
        _logger.info('Runing sextractor in %s', basename)
        sex.run('%s,%s' % (basename, basename))
        
        # Sort catalog by flux
        catalog = sex.catalog()
        catalog = sorted(catalog, key=operator.itemgetter('FLUX_BEST'), reverse=True)
        
        # set of indices of the N first objects
        OBJS_I_KEEP = 3
        indices = set(obj['NUMBER'] for obj in catalog[:OBJS_I_KEEP])
        
        base = numpy.empty((len(images_info), OBJS_I_KEEP))
        error = numpy.empty((len(images_info), OBJS_I_KEEP))
        
        for idx, image in enumerate(images_info):
            imagename = name_skysub_proc(image.baselabel, step)

            sex.config['CATALOG_NAME'] = 'catalogue-%s-i%01d.cat' % (image.baselabel, step)

            # Lauch SExtractor on a FITS file
            # om double image mode
            _logger.info('Runing sextractor in %s', imagename)
            sex.run('%s,%s' % (basename, imagename))
            catalog = sex.catalog()
            
            # Extinction correction
            excor = pow(10, -0.4 * image.airmass * self.parameters['extinction'])
            excor = 1.0
            base[idx] = [obj['FLUX_BEST'] / excor
                                     for obj in catalog if obj['NUMBER'] in indices]
            error[idx] = [obj['FLUXERR_BEST'] / excor
                                     for obj in catalog if obj['NUMBER'] in indices]
        
        data = base / base[0]
        err = error / base[0] # sigma
        w = 1 / err / err
        # weighted mean of the flux values
        wdata = numpy.average(data, axis=1, weights=w)
        wsigma = 1 / numpy.sqrt(w.sum(axis=1))
        
        # Actions to carry over images when checking the flux
        # of the objects in different images
        def warn_action(img):
            _logger.warn('Image %s has low flux in objects', img.baselabel)
            img.valid_science = True
        
        def reject_action(img):
            img.valid_science = False
            _logger.info('Image %s rejected, has low flux in objects', img.baselabel)            
            pass
        
        def default_action(img):
            _logger.info('Image %s accepted, has correct flux in objects', img.baselabel)      
            img.valid_science = True
        
        # Actions
        dactions = {'warn': warn_action, 'reject': reject_action, 'default': default_action}

        levels = self.parameters['check_photometry_levels']
        actions = self.parameters['check_photometry_actions']
        
        x = range(len(images_info))
        vals, (_, sigma) = self.check_photometry_categorize(x, wdata, 
                                                           levels, tags=actions)
        # n sigma level to plt
        n = 3
        self.check_photometry_plot(vals, wsigma, levels, n * sigma)
        
        for x, _, t in vals:
            try:
                action = dactions[t]
            except KeyError:
                _logger.warning('Action named %s not recognized, ignoring', t)
                action = default_action
            for p in x:
                action(images_info[p])                
                
    def check_position(self, images_info, sf_data, seeing_fwhm):

        _logger.info('Checking positions')
        # Check position of bright objects
        weigthmap = 'weights4rms.fits'
        
        wmap = numpy.zeros_like(sf_data[0])
        
        # Center of the image
        border = 300
        wmap[border:-border, border:-border] = 1                    
        pyfits.writeto(weigthmap, wmap, clobber=True)
        
        basename = 'result_i%0d.fits' % (step)
        sex = SExtractor()
        sex.config['VERBOSE_TYPE'] = 'QUIET'
        sex.config['PIXEL_SCALE'] = 1
        sex.config['BACK_TYPE'] = 'AUTO'
        if  seeing_fwhm is not None and seeing_fwhm > 0:
            sex.config['SEEING_FWHM'] = seeing_fwhm * sex.config['PIXEL_SCALE']
        sex.config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
        sex.config['WEIGHT_IMAGE'] = weigthmap
        
        sex.config['PARAMETERS_LIST'].append('FLUX_BEST')
        sex.config['PARAMETERS_LIST'].append('FLUXERR_BEST')
        sex.config['PARAMETERS_LIST'].append('FWHM_IMAGE')
        sex.config['PARAMETERS_LIST'].append('CLASS_STAR')
        
        sex.config['CATALOG_NAME'] = 'master-catalogue-i%01d.cat' % step
        
        _logger.info('Runing sextractor in %s', basename)
        sex.run('%s,%s' % (basename, basename))
        
        # Sort catalog by flux
        catalog = sex.catalog()
        catalog = sorted(catalog, key=operator.itemgetter('FLUX_BEST'), reverse=True)
        
        # set of indices of the N first objects
        OBJS_I_KEEP = 10
        
        master = [(obj['X_IMAGE'], obj['Y_IMAGE']) for obj in catalog[:OBJS_I_KEEP]]
        
        for image in images_info:
            imagename = name_skysub_proc(image.baselabel, self.iter)

            sex.config['CATALOG_NAME'] = 'catalogue-self-%s-i%01d.cat' % (image.baselabel, step)

            # Lauch SExtractor on a FITS file
            # on double image mode
            _logger.info('Runing sextractor in %s', imagename)
            sex.run(imagename)
            catalog = sex.catalog()
            
            
            data = [(obj['X_IMAGE'], obj['Y_IMAGE']) for obj in catalog]
            
            tree = KDTree(data)
            
            # Search 2 neighbors
            dists, _ids = tree.query(master, 2, distance_upper_bound=5)
            
            for i in dists[:,0]:
                print i
            
            
            _logger.info('Mean offset correction for image %s is %f', imagename, dists[:,0].mean())
            #raw_input('press any key')
                            
                
                
    def check_photometry_categorize(self, x, y, levels, tags=None):
        '''Put every point in its category.
    
        levels must be sorted.'''   
        x = numpy.asarray(x)
        y = numpy.asarray(y)
        ys = y.copy()
        ys.sort()
        # Mean of the upper half
        m = ys[len(ys) / 2:].mean()
        y /= m
        m = 1.0
        s = ys[len(ys) / 2:].std()
        result = []

        if tags is None:
            tags = range(len(levels) + 1)

        for l, t in zip(levels, tags):
            indc = y < l
            if indc.any():
                x1 = x[indc]
                y1 = y[indc]
                result.append((x1, y1, t))

                x = x[indc == False]
                y = y[indc == False]
        else:
            result.append((x, y, tags[-1]))

        return result, (m,s)     

    def create_mask(self, sf_data, seeing_fwhm):
        # FIXME more plots
        self.figure_final_before_s(sf_data[0])

        #
        remove_border = True
        
        # sextractor takes care of bad pixels
        sex = SExtractor()
        sex.config['CHECKIMAGE_TYPE'] = "SEGMENTATION"
        sex.config["CHECKIMAGE_NAME"] = name_segmask(step)
        sex.config['VERBOSE_TYPE'] = 'QUIET'
        sex.config['PIXEL_SCALE'] = 1
        sex.config['BACK_TYPE'] = 'AUTO' 

        if seeing_fwhm is not None and seeing_fwhm > 0:
            sex.config['SEEING_FWHM'] = seeing_fwhm * sex.config['PIXEL_SCALE']

        sex.config['PARAMETERS_LIST'].append('FLUX_BEST')
        sex.config['PARAMETERS_LIST'].append('X_IMAGE')
        sex.config['PARAMETERS_LIST'].append('Y_IMAGE')
        sex.config['PARAMETERS_LIST'].append('A_IMAGE')
        sex.config['PARAMETERS_LIST'].append('B_IMAGE')
        sex.config['PARAMETERS_LIST'].append('THETA_IMAGE')
        sex.config['PARAMETERS_LIST'].append('FWHM_IMAGE')
        sex.config['PARAMETERS_LIST'].append('CLASS_STAR')                
        if remove_border:
            weigthmap = 'weights4rms.fits'
            
            # Create weight map, remove n pixs from either side
            # using a Hannig filter
            # npix = 90          
            # w1 = npix
            # w2 = npix
            # wmap = numpy.ones_like(sf_data[0])
            
            # cos_win1 = numpy.hanning(2 * w1)
            # cos_win2 = numpy.hanning(2 * w2)
                                   
            # wmap[:,:w1] *= cos_win1[:w1]                    
            # wmap[:,-w1:] *= cos_win1[-w1:]
            # wmap[:w2,:] *= cos_win2[:w2, numpy.newaxis]
            # wmap[-w2:,:] *= cos_win2[-w2:, numpy.newaxis]                 
            
            # Take the number of combined images from the combined image
            wm = sf_data[2].copy()
            # Dont search objects where nimages < lower
            # FIXME: this is a magic number
            # We ignore objects in regions where we have less
            # than 10% of the images
            lower = sf_data[2].max() / 10
            wm[wm < lower] = 0
            pyfits.writeto(weigthmap, wm, clobber=True)
                                
            sex.config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            # FIXME: this is a magic number
            # sex.config['WEIGHT_THRESH'] = 50
            sex.config['WEIGHT_IMAGE'] = weigthmap
        
        filename = 'result_i%0d.fits' % (step)
        
        # Lauch SExtractor on a FITS file
        sex.run(filename)
        
        # Plot objects
        # FIXME, plot sextractor objects on top of image
        patches = []
        fwhms = []
        nfirst = 0
        catalog_f = sopen(sex.config['CATALOG_NAME'])
        try:
            star = catalog_f.readline()
            while star:
                flags = star['FLAGS']
                # ignoring those objects with corrupted apertures
                if flags & sexcatalog.CORRUPTED_APER:
                    star = catalog_f.readline()
                    continue
                center = (star['X_IMAGE'], star['Y_IMAGE'])
                wd = 10 * star['A_IMAGE']
                hd = 10 * star['B_IMAGE']
                color = 'red'
                e = Ellipse(center, wd, hd, star['THETA_IMAGE'], color=color)
                patches.append(e)
                fwhms.append(star['FWHM_IMAGE'])
                nfirst += 1
                # FIXME Plot a ellipse
                star = catalog_f.readline()
        finally:
            catalog_f.close()
            
        p = PatchCollection(patches, alpha=0.4)
        ax = self._figure.gca()
        ax.add_collection(p)
        self._figure.canvas.draw()
        self._figure.savefig('figure-segmentation-overlay_%01d.png' % step)

        self.figure_fwhm_histogram(fwhms)
                    
        # mode with an histogram
        hist, edges = numpy.histogram(fwhms, 50)
        idx = hist.argmax()
        
        seeing_fwhm = 0.5 * (edges[idx] + edges[idx + 1]) 
        if seeing_fwhm <= 0:
            _logger.warning('Seeing FHWM %f pixels is negative, reseting', seeing_fwhm)
            seeing_fwhm = None
        else:
            _logger.info('Seeing FHWM %f pixels (%f arcseconds)', seeing_fwhm, seeing_fwhm * sex.config['PIXEL_SCALE'])
        objmask = pyfits.getdata(name_segmask(step))
        return objmask, seeing_fwhm
    

    def resize(self, images_info):
        _logger.info('Computing offsets')
        image_shapes = images_info[0].baseshape
        offsets = [image.offset for image in images_info if image.valid_science]        
        finalshape, offsetsp = combine_shape(image_shapes, offsets)
        _logger.info('Shape of resized array is %s', finalshape)

        _logger.info('Resizing images and masks')            
        
        for image, noffset in zip(images_info, offsetsp):
            if image.valid_science:
                region, _ = subarray_match(finalshape, noffset, image.baseshape)
            else:
                region = None
                noffset = None
            image.valid_region = region
            image.noffset = noffset
            imgn, maskn = name_redimensioned_images(image.baselabel, step)
            image.resized_base = imgn
            image.resized_mask = maskn
                    
            self.resize_frame_and_mask(image, finalshape, imgn, maskn, scale=1)

        return images_info
            
    def run(self, obresult):
        
        metadata = self.instrument['metadata']
        detector = self.instrument['detectors'][0] # One detector
        amplifiers = self.instrument['amplifiers'][0] # Amplifiers
        
        amplifiers = [(slice(chan[0][0], chan[0][1]), slice(chan[1][0], chan[1][1]))  
                      for chan in amplifiers]
        
        # States
        BASIC, PRERED, CHECKRED, FULLRED, COMPLETE = range(5)

        state = BASIC
        
        images_info = []
        
        # Creating ImageInformation
        for image, offset in obresult.frames:
            # Label
            ii = ImageInformation()
            ii.baselabel = os.path.splitext(image)[0]
            ii.label = image
            ii.base = image
            ii.mask = self.parameters['master_bpm']
            ii.offset = offset
 
            ii.valid_region = None
            ii.resized_base = None
            ii.resized_mask = None
            ii.objmask = None
            ii.objmask_data = None
            hdr = pyfits.getheader(ii.base)
            try:
                ii.exposure = hdr[str(metadata['exposure'])]
                ii.baseshape = get_image_shape(hdr)
                ii.airmass = hdr[str(metadata['airmass'])]
                ii.mjd = hdr[str(metadata['juliandate'])]
            except KeyError, e:
                raise RecipeError("%s in image %s" % (str(e), ii.base))
            images_info.append(ii)
    
        images_info = update_sky_related(images_info, nimages=self.parameters['sky_images'])
        image_shapes = images_info[0].baseshape
    
        niteration = self.parameters['iterations']
        # Final image, not yet built
        sf_data = None
        seeing_fwhm = None
        self.iter = 0

        while True:
            if state == BASIC:    
                _logger.info('Basic processing')

                # Basic processing
                # Open bpm, dark and flat
                bpm = pyfits.getdata(self.parameters['master_bpm'])
                mdark = pyfits.getdata(self.parameters['master_dark'])
                # FIXME: Use a flat field image eventually
                #mflat = pyfits.getdata(self.parameters['master_flat'].filename)
        
                basicflow = SerialFlow([BadPixelCorrector(bpm),
                                NonLinearityCorrector(self.parameters['nonlinearity']),
                                DarkCorrector(mdark)])

                for image in images_info:
                    self.basic_processing(image, basicflow)
        
                del bpm
                del mdark
                del basicflow
                state = PRERED

            elif state == PRERED:
                # Resizing images                            
                images_info = self.resize(images_info)

                # startup plot system



                # superflat
                _logger.info('Iter %d, superflat correction (SF)', self.iter)
                # Compute scale factors (median)           
                images_info = self.update_scale_factors(images_info, self.iter)

                # Create superflat
                superflat = self.compute_superflat(images_info, None, amplifiers)
            
                # Apply superflat
                self.figure_init(images_info[0].baseshape)
                images_info = self.apply_superflat(images_info, superflat, step=self.iter, save=True)

                _logger.info('Simple sky correction')
                for image in images_info:            
                    self.compute_simple_sky(image, self.iter, save=True)
                
                # Combining the images
                _logger.info("Iter %d, Combining the images", self.iter)
                # FIXME: only for science
                sf_data = self.combine_frames(images_info, self.iter)
            
                self.figures_after_combine(sf_data)
                      
                _logger.info('Iter %d, finished', self.iter)

                state = CHECKRED
            elif state == CHECKRED:
                
                seeing_fwhm = None

                #self.check_position(images_info, sf_data, seeing_fwhm)
                recompute = False
                if recompute:
                    _logger.info('Recentering is needed')
                    state = PRERED
                else:
                    _logger.info('Recentering is not needed')
                    _logger.info('Checking photometry')
                    self.check_photometry(images_info, sf_data, seeing_fwhm)
                    state = FULLRED
                    
            elif state == FULLRED:

                # Generating segmentation image
                _logger.info('Iter %d, generating segmentation image', self.iter)
                objmask, seeing_fwhm = self.create_mask(sf_data, seeing_fwhm)
                self.iter +=1
                # Update objects mask
                # For all images    
                for image in images_info:
                    image.objmask = name_object_mask(image.baselabel, self.iter)
                    _logger.info('Iter %d, create object mask %s', self.iter,  image.objmask)                 
                    image.objmask_data = objmask[image.valid_region]
                    pyfits.writeto(image.objmask, image.objmask_data, clobber=True)

                _logger.info('Iter %d, superflat correction (SF)', self.iter)
                # Compute scale factors (median)           
                images_info = self.update_scale_factors(images_info, self.iter)
            
                # Combining images to obtain the sky flat
                superflat = self.compute_superflat(images_info, objmask, amplifiers)
    
                # Apply superflat
                self.figure_init(images_info[0].baseshape)
                
                images_info = self.apply_superflat(images_info, superflat, step=self.iter, save=True)

                _logger.info('Iter %d, advanced sky correction (SC)', self.iter)
                # FIXME: Only for science          
                for image in images_info:
                    if image.valid_science:       
                        self.compute_advanced_sky(image, objmask)
            
                # Combining the images
                _logger.info("Iter %d, Combining the images", self.iter)
                # FIXME: only for science
                sf_data = self.combine_frames(images_info, self.iter)
                self.figures_after_combine(sf_data)

                if self.iter >= niteration:
                    state = COMPLETE
            else:
                break

        primary_headers = {'FILENAME': self.parameters['resultname'],
                           }
        result = create_result(sf_data[0], headers=primary_headers,
                                variance=sf_data[1], 
                                exmap=sf_data[2].astype('int16'))
        
        _logger.info("Final image created")
        return {'products': [DataFrame(result), SourcesCatalog()]}

