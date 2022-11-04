#
# Copyright 2011-2017 Universidad Complutense de Madrid
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

"""Routines shared by image mode recipes"""

import os
import logging
import shutil
import math

import six
import numpy
from astropy.io import fits
import astropy.wcs
from astropy.visualization import SqrtStretch
from astropy.visualization import PercentileInterval
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy.spatial import KDTree as KDTree
import sep

import matplotlib as mpl
from matplotlib.patches import Ellipse
from matplotlib.collections import PatchCollection

from numina import __version__
from numina.core import DataFrame
from numina.array import fixpix2
from numina.frame import resize_fits, custom_region_to_str
from numina.array import combine_shape, correct_flatfield
from numina.array import subarray_match
from numina.array.combine import flatcombine, median, quantileclip
from numina.exceptions import RecipeError

from emirdrp.util.sextractor import SExtractor
from emirdrp.util.sextractor import open as sopen
from emirdrp.core.recipe import EmirRecipe
from emirdrp.products import SourcesCatalog
from emirdrp.instrument.channels import FULL
from emirdrp.processing.wcs import offsets_from_wcs

from .checks import check_photometry
from .naming import (name_redimensioned_frames, name_object_mask,
                     name_skybackground)
from .naming import name_skybackgroundmask, name_skysub_proc, name_skyflat
from .naming import name_skyflat_proc, name_segmask


_logger = logging.getLogger('numina.recipes.emir')


def intersection(a, b, scale=1):
    '''Intersection between two segments.'''
    try:
        a1, a2 = a
    except TypeError:
        a1 = a.start
        a2 = a.stop

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


class DirectImageCommon(EmirRecipe):

    logger = _logger
    BASIC, PRERED, CHECKRED, FULLRED, COMPLETE = [0, 1, 2, 3, 4]
    __version__ = '1'

    def __init__(self, *args, **kwds):
        super(DirectImageCommon, self).__init__(version=__version__)
        # Required to delay the backend initialization (issue numina/#102, #49)
        import matplotlib.pyplot as plt

        self._figure = plt.figure(facecolor='white')
        # FigureCanvasAgg has no set_window_title
        # self._figure.canvas.set_window_title('Recipe Plots')
        self._figure.canvas.draw()

    def process(self, ri,
                window=None, subpix=1,
                store_intermediate=True,
                target_is_sky=True, stop_after=PRERED):

        numpy.seterr(divide='raise')

        recipe_input = ri
        # FIXME: hardcoded instrument information
        keywords = {'airmass': 'AIRMASS',
                    'exposure': 'EXPTIME',
                    'imagetype': 'IMGTYP',
                    'juliandate': 'MJD-OBS',
                    'tstamp': 'TSTAMP'
                    }
        baseshape = [2048, 2048]
        channels = FULL

        if window is None:
            window = tuple((0, siz) for siz in baseshape)

        if store_intermediate:
            pass

        # States
        sf_data = None
        state = self.BASIC
        step = 0

        try:
            niteration = ri.iterations
        except KeyError:
            niteration = 1

        while True:
            if state == self.BASIC:
                _logger.info('Basic processing')

                # Basic processing
                basicflow = self.init_filters(recipe_input)

                for frame in ri.obresult.frames:
                    with frame.open() as hdulist:
                        hdulist = basicflow(hdulist)

                if stop_after == state:
                    break
                else:
                    state = self.PRERED
            elif state == self.PRERED:
                # Shape of the window
                windowshape = tuple((i[1] - i[0]) for i in window)
                _logger.debug('Shape of window is %s', windowshape)
                # Shape of the scaled window
                subpixshape = tuple((side * subpix) for side in windowshape)

                # Scaled window region
                scalewindow = tuple(
                    slice(*(subpix * i for i in p)) for p in window)
                # Window region
                window = tuple(slice(*p) for p in window)

                scaled_chan = clip_slices(channels, window, scale=subpix)

                # Reference pixel in the center of the frame
                refpix = numpy.divide(
                    numpy.array([baseshape], dtype='int'), 2).astype('float')

                # lists of targets and sky frames
                targetframes = []
                skyframes = []

                for frame in ri.obresult.frames:

                    # Getting some metadata from FITS header
                    hdr = fits.getheader(frame.label)
                    try:
                        frame.exposure = hdr[str(keywords['exposure'])]
                        # frame.baseshape = get_image_shape(hdr)
                        frame.airmass = hdr[str(keywords['airmass'])]
                        frame.mjd = hdr[str(keywords['tstamp'])]
                    except KeyError as e:
                        raise KeyError("%s in frame %s" %
                                       (str(e), frame.label))

                    frame.baselabel = os.path.splitext(frame.label)[0]
                    frame.mask = ri.master_bpm
                    # Insert pixel offsets between frames
                    frame.objmask_data = None
                    frame.valid_target = False
                    frame.valid_sky = False
                    frame.valid_region = scalewindow
                    # FIXME: hardcode itype for the moment
                    frame.itype = 'TARGET'
                    if frame.itype == 'TARGET':
                        frame.valid_target = True
                        targetframes.append(frame)
                        if target_is_sky:
                            frame.valid_sky = True
                            skyframes.append(frame)
                    if frame.itype == 'SKY':
                        frame.valid_sky = True
                        skyframes.append(frame)

#                labels = [frame.label for frame in targetframes]

                if ri.offsets is not None:
                    _logger.info('Using offsets from parameters')
                    base_ref = numpy.asarray(ri.offsets)
                    list_of_offsets= -(base_ref - base_ref[0])
                else:
                    _logger.info('Computing offsets from WCS information')
                    list_of_offsets = offsets_from_wcs(targetframes, refpix)

                # FIXME: Im using offsets in row/columns
                # the values are provided in XY so flip-lr
                list_of_offsets = numpy.fliplr(list_of_offsets)

                # Insert pixel offsets between frames
                for frame, off in zip(targetframes, list_of_offsets):

                    # Insert pixel offsets between frames
                    frame.pix_offset = off
                    frame.scaled_pix_offset = subpix * off

                    _logger.debug('Frame %s, offset=%s, scaled=%s',
                                  frame.label, off, subpix * off)

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
                self.update_scale_factors(ri.obresult.frames)

                # Create superflat
                superflat = self.compute_superflat(skyframes,
                                                   channels=scaled_chan,
                                                   step=step)

                # Apply superflat
                self.figure_init(subpixshape)
                self.apply_superflat(ri.obresult.frames, superflat)

                _logger.info('Simple sky correction')
                if target_is_sky:
                    # Each frame is the closest sky frame available

                    for frame in ri.obresult.frames:
                        self.compute_simple_sky_for_frame(frame, frame)
                else:
                    self.compute_simple_sky(targetframes, skyframes)

                # Combining the frames
                _logger.info("Step %d, Combining target frames", step)

                sf_data = self.combine_frames(
                    targetframes, extinction=ri.extinction)

                self.figures_after_combine(sf_data)

                _logger.info('Step %d, finished', step)

                if stop_after == state:
                    break
                else:
                    state = self.CHECKRED
            elif state == self.CHECKRED:

                seeing_fwhm = None

                # self.check_position(images_info, sf_data, seeing_fwhm)
                recompute = False
                if recompute:
                    _logger.info('Recentering is needed')
                    state = self.PRERED
                else:
                    _logger.info('Recentering is not needed')
                    _logger.info('Checking photometry')
                    check_photometry(targetframes, sf_data,
                                     seeing_fwhm, figure=self._figure)

                    if stop_after == state:
                        break
                    else:
                        state = self.FULLRED
            elif state == self.FULLRED:

                # Generating segmentation image
                _logger.info('Step %d, generating segmentation image', step)
                objmask, seeing_fwhm = self.create_mask(
                    sf_data, seeing_fwhm, step=step)
                step += 1
                # Update objects mask
                # For all images
                # FIXME:
                for frame in targetframes:
                    frame.objmask = name_object_mask(frame.baselabel, step)
                    _logger.info(
                        'Step %d, create object mask %s', step,  frame.objmask)
                    frame.objmask_data = objmask[frame.valid_region]
                    fits.writeto(
                        frame.objmask, frame.objmask_data, overwrite=True)

                if not target_is_sky:
                    # Empty object mask for sky frames
                    bogus_objmask = numpy.zeros(windowshape, dtype='int')

                    for frame in skyframes:
                        frame.objmask_data = bogus_objmask

                _logger.info('Step %d, superflat correction (SF)', step)

                # Compute scale factors (median)
                self.update_scale_factors(ri.obresult.frames, step)

                # Create superflat
                superflat = self.compute_superflat(skyframes, scaled_chan,
                                                   segmask=objmask, step=step)

                # Apply superflat
                self.figure_init(subpixshape)

                self.apply_superflat(
                    ri.obresult.frames, superflat, step=step, save=True)

                _logger.info('Step %d, advanced sky correction (SC)', step)
                self.compute_advanced_sky(targetframes, objmask,
                                          skyframes=skyframes,
                                          target_is_sky=target_is_sky,
                                          step=step)

                # Combining the images
                _logger.info("Step %d, Combining the images", step)
                # FIXME: only for science
                sf_data = self.combine_frames(
                    targetframes, ri.extinction, step=step)
                self.figures_after_combine(sf_data)

                if step >= niteration:
                    state = self.COMPLETE
            else:
                break

        if sf_data is None:
            raise RecipeError(
                'no combined image has been generated at step %d', state)

        hdu = fits.PrimaryHDU(sf_data[0])
        hdr = hdu.header
        hdr.update('NUMXVER', __version__, 'Numina package version')
        hdr.update('NUMRNAM', self.__class__.__name__, 'Numina recipe name')
        hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

        hdr.update('FILENAME', 'result.fits')
        hdr.update('IMGTYP', 'TARGET', 'Image type')
        hdr.update('NUMTYP', 'TARGET', 'Data product type')

        varhdu = fits.ImageHDU(sf_data[1], name='VARIANCE')
        num = fits.ImageHDU(sf_data[2], name='MAP')

        result = fits.HDUList([hdu, varhdu, num])

        _logger.info("Final frame created")

        return DataFrame(result), SourcesCatalog()

    def compute_simple_sky(self, targetframes, skyframes,
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
        _dis, idxs = kdtree.query(tarray, k=1,
                                  distance_upper_bound=maxsep * MIN_TO_DAY)

        for tid, idss in enumerate(idxs):
            try:
                tf = targetframes[tid]
                sf = skyframes[idss]
                self.compute_simple_sky_for_frame(tf, sf, step=step, save=save)
            except IndexError:
                _logger.error(
                    'No sky image available for frame %s', tf.lastname)
                raise

    def compute_simple_sky_for_frame(self, frame, skyframe, step=0, save=True):
        _logger.info('Correcting sky in frame %s', frame.lastname)
        _logger.info('with sky computed from frame %s', skyframe.lastname)

        if hasattr(skyframe, 'median_sky'):
            sky = skyframe.median_sky
        else:

            with fits.open(skyframe.lastname, mode='readonly') as hdulist:
                data = hdulist['primary'].data
                valid = data[frame.valid_region]

                if skyframe.objmask_data is not None:
                    _logger.debug('object mask defined')
                    msk = frame.objmask_data
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

        with fits.open(frame.lastname, mode='update') as hdulist:
            data = hdulist['primary'].data
            valid = data[frame.valid_region]
            valid -= sky

    def compute_advanced_sky(self, targetframes, objmask,
                             skyframes=None, target_is_sky=False,
                             maxsep=5.0,
                             nframes=10,
                             step=0, save=True):

        if target_is_sky:
            skyframes = targetframes
            # Each frame is its closets sky frame
            nframes += 1
        elif skyframes is None:
            raise ValueError('skyframes not defined')

        # build kdtree
        sarray = numpy.array([frame.mjd for frame in skyframes])
        # shape must be (n, 1)
        sarray = numpy.expand_dims(sarray, axis=1)

        # query
        tarray = numpy.array([frame.mjd for frame in targetframes])
        # shape must be (n, 1)
        tarray = numpy.expand_dims(tarray, axis=1)

        kdtree = KDTree(sarray)

        # 1 / minutes in a Julian day
        SCALE = 60.0
        # max_time_sep = ri.sky_images_sep_time / 1440.0
        _dis, idxs = kdtree.query(tarray, k=nframes,
                                  distance_upper_bound=maxsep * SCALE)

        nsky = len(sarray)

        for tid, idss in enumerate(idxs):
            try:
                tf = targetframes[tid]
                _logger.info('Step %d, SC: computing advanced sky for %s',
                             step, tf.baselabel)
                # filter(lambda x: x < nsky, idss)
                locskyframes = []
                for si in idss:
                    if tid == si:
                        # this sky frame its the current frame, reject
                        continue
                    if si < nsky:
                        _logger.debug('Step %d, SC: %s is a sky frame',
                                      step, skyframes[si].baselabel)
                        locskyframes.append(skyframes[si])
                self.compute_advanced_sky_for_frame(
                    tf, locskyframes, step=step, save=save)
            except IndexError:
                _logger.error(
                    'No sky image available for frame %s', tf.lastname)
                raise

    def compute_advanced_sky_for_frame(self, frame, skyframes,
                                       step=0, save=True):
        _logger.info('Correcting sky in frame %s', frame.lastname)
        _logger.info('with sky computed from frames')
        for i in skyframes:
            _logger.info('%s', i.flat_corrected)

        data = []
        scales = []
        masks = []
        # handle the FITS file to close it finally
        desc = []
        try:
            for i in skyframes:
                filename = i.flat_corrected
                hdulist = fits.open(filename, mode='readonly', memmap=True)

                data.append(hdulist['primary'].data[i.valid_region])
                desc.append(hdulist)
                scales.append(numpy.median(data[-1]))
                if i.objmask_data is not None:
                    masks.append(i.objmask_data)
                    _logger.debug('object mask is shared')
                elif i.objmask is not None:
                    hdulistmask = fits.open(
                        i.objmask, mode='readonly', memmap=True)
                    masks.append(hdulistmask['primary'].data)
                    desc.append(hdulistmask)
                    _logger.debug('object mask is particular')
                else:
                    _logger.warn('no object mask for %s', filename)

            _logger.debug('computing background with %d frames', len(data))
            sky, _, num = median(data, masks, scales=scales)

        finally:
            # Closing all FITS files
            for hdl in desc:
                hdl.close()

        if numpy.any(num == 0):
            # We have pixels without
            # sky background information
            _logger.warn('pixels without sky information when correcting %s',
                         frame.flat_corrected)
            binmask = num == 0
            # FIXME: during development, this is faster
            # sky[binmask] = sky[num != 0].mean()

            # To continue we interpolate over the patches
            fixpix2(sky, binmask, out=sky, iterations=1)
            name = name_skybackground(frame.baselabel, step)
            fits.writeto(name, sky, overwrite=True)
            name = name_skybackgroundmask(frame.baselabel, step)
            fits.writeto(name, binmask.astype('int16'), overwrite=True)

        dst = name_skysub_proc(frame.baselabel, step)
        prev = frame.lastname
        shutil.copyfile(prev, dst)
        frame.lastname = dst

        with fits.open(frame.lastname, mode='update') as hdulist:
            data = hdulist['primary'].data
            valid = data[frame.valid_region]
            valid -= sky

    def combine_frames(self, frames, extinction, out=None, step=0):
        _logger.debug('Step %d, opening sky-subtracted frames', step)

        def fits_open(name):
            '''Open FITS with memmap in readonly mode'''
            return fits.open(name, mode='readonly', memmap=True)

        frameslll = [fits_open(frame.lastname)
                     for frame in frames if frame.valid_target]
        _logger.debug('Step %d, opening mask frames', step)
        mskslll = [fits_open(frame.resized_mask)
                   for frame in frames if frame.valid_target]
        _logger.debug('Step %d, combining %d frames', step, len(frameslll))
        try:
            extinc = [pow(10, -0.4 * frame.airmass * extinction)
                      for frame in frames if frame.valid_target]
            data = [i['primary'].data for i in frameslll]
            masks = [i['primary'].data for i in mskslll]

            out = quantileclip(data, masks, scales=extinc,
                               dtype='float32', out=out, fclip=0.1)

            # saving the three extensions
            fits.writeto('result_i%0d.fits' % step, out[0], overwrite=True)
            fits.writeto('result_var_i%0d.fits' % step, out[1], overwrite=True)
            fits.writeto('result_npix_i%0d.fits' % step, out[2], overwrite=True)

            return out

        finally:
            _logger.debug('Step %d, closing sky-subtracted frames', step)
            for f in frameslll:
                f.close()
            _logger.debug('Step %d, closing mask frames', step)
            for f in mskslll:
                f.close()

    def apply_superflat(self, frames, flatdata, step=0, save=True):
        _logger.info("Step %d, SF: apply superflat", step)

        # Process all frames with the fitted flat
        # FIXME: not sure
        for frame in frames:
            self.correct_superflat(frame, flatdata, step=step, save=save)
        return frames

    def correct_superflat(self, frame, fitted, step=0, save=True):

        frame.flat_corrected = name_skyflat_proc(frame.baselabel, step)

        if save:
            shutil.copyfile(frame.resized_base, frame.flat_corrected)
        else:
            os.rename(frame.resized_base, frame.flat_corrected)

        _logger.info("Step %d, SF: apply superflat to frame %s",
                     step, frame.flat_corrected)
        with fits.open(frame.flat_corrected, mode='update') as hdulist:
            data = hdulist['primary'].data
            datar = data[frame.valid_region]
            data[frame.valid_region] = correct_flatfield(datar, fitted)

            frame.lastname = frame.flat_corrected

            # FIXME: plotting
            try:
                self.figure_image(data[frame.valid_region], frame)
            except ValueError:
                _logger.warning('Problem plotting %s', frame.lastname)

    def compute_superflat(self, frames, channels, segmask=None, step=0):
        _logger.info("Step %d, SF: combining the frames without offsets", step)
        try:
            filelist = []
            data = []
            masks = []
            for frame in frames:
                _logger.debug('Step %d, opening resized frame %s',
                              step, frame.resized_base)
                hdulist = fits.open(
                    frame.resized_base, memmap=True, mode='readonly')
                filelist.append(hdulist)
                data.append(hdulist['primary'].data[frame.valid_region])

            scales = [frame.median_scale for frame in frames]

            # FIXME: plotting
            self.figure_median_background(scales)

            if segmask is not None:
                masks = [segmask[frame.valid_region] for frame in frames]
            else:
                for frame in frames:
                    _logger.debug('Step %d, opening resized mask %s',
                                  step, frame.resized_mask)
                    hdulist = fits.open(
                        frame.resized_mask, memmap=True, mode='readonly')
                    filelist.append(hdulist)
                    masks.append(hdulist['primary'].data[frame.valid_region])

            _logger.debug('Step %d, combining %d frames', step, len(data))

            sf_data, _sf_var, sf_num = flatcombine(data, masks, scales=scales,
                                                   blank=1.0 / scales[0])
        finally:
            _logger.debug('Step %d, closing resized frames and mask', step)
            for fileh in filelist:
                fileh.close()

        # We interpolate holes by channel
        _logger.debug('Step %d, interpolating holes by channel', step)
        for channel in channels:
            mask = (sf_num[channel] == 0)
            if numpy.any(mask):
                fixpix2(sf_data[channel], mask, out=sf_data[channel])

        # Normalize, flat has mean = 1
        sf_data /= sf_data.mean()

        # Auxiliary data
        sfhdu = fits.PrimaryHDU(sf_data)
        sfhdu.writeto(name_skyflat('comb', step), overwrite=True)
        return sf_data

    def update_scale_factors(self, frames, step=0):
        _logger.info('Step %d, SF: computing scale factors', step)
        # FIXME: not sure
        for frame in frames:
            region = frame.valid_region
            data = fits.getdata(frame.resized_base)[region]
            mask = fits.getdata(frame.resized_mask)[region]
            # FIXME: while developing this ::10 is faster, remove later
            frame.median_scale = numpy.median(data[mask == 0][::10])
            _logger.debug('median value of %s is %f',
                          frame.resized_base, frame.median_scale)
        return frames

    def resize(self, frames, shape, offsetsp, finalshape, window=None,
               scale=1, step=0):
        _logger.info('Resizing frames and masks')
        for frame, rel_offset in zip(frames, offsetsp):
            if frame.valid_target:
                region, _ = subarray_match(finalshape, rel_offset, shape)
                # Valid region
                frame.valid_region = region
                # Relative offset
                frame.rel_offset = rel_offset
                # names of frame and mask
                framen, maskn = name_redimensioned_frames(
                    frame.baselabel, step)
                frame.resized_base = framen
                frame.resized_mask = maskn
                _logger.debug('%s, valid region is %s, relative offset is %s',
                              frame.label, custom_region_to_str(region),
                              rel_offset
                              )
                self.resize_frame_and_mask(
                    frame, finalshape, framen, maskn, window, scale)

        return frames

    def resize_frame_and_mask(self, frame, finalshape,
                              framen, maskn, window, scale):
        _logger.info('Resizing frame %s, window=%s, subpix=%i', frame.label,
                     custom_region_to_str(window), scale)
        hdul = frame.open()
        baseshape = hdul[0].data.shape

        # FIXME: Resize_fits saves the resized image in framen
        resize_fits(hdul, framen, finalshape, frame.valid_region,
                    window=window, scale=scale, dtype='float32')

        _logger.info('Resizing mask %s, subpix x%i', frame.label, scale)
        # We don't conserve the sum of the values of the frame here, just
        # expand the mask
        if frame.mask is None:
            self.logger.warning('BPM missing, use zeros instead')
            false_mask = numpy.zeros(baseshape, dtype='int16')
            hdum = fits.HDUList(fits.PrimaryHDU(false_mask))
            frame.mask = DataFrame(frame=hdum)
        resize_fits(frame.mask.open(), maskn, finalshape, frame.valid_region,
                    fill=1, window=window, scale=scale, conserve=False)

    def figure_init(self, shape):
        self._figure.clf()
        ax = self._figure.add_subplot(111)
        cmap = mpl.cm.get_cmap('gray')
        norm = mpl.colors.LogNorm()
        ax.imshow(numpy.ones(shape), cmap=cmap, norm=norm)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        self._figure.canvas.draw()

    def figures_after_combine(self, sf_data, step=0):

        # FIXME, more plots
        def truncated(array, frac=0.1):
            '''Dirty truncated mean'''
            nf = int(array.size * frac)
            array.sort()
            new = array[nf:-nf]
            lnew = len(new)
            if lnew == 0:
                return 0, 0
            elif lnew == 1:
                return new.mean(), 0.0
            else:
                return new.mean(), new.std()

        ndata = sf_data[2].astype('int')
        data = sf_data[0]

        nimages = ndata.max()

        rnimage = list(six.moves.range(1, nimages + 1))
        rmean = rnimage[:]
        rstd = rnimage[:]

        for pix in rnimage:
            rmean[pix - 1], rstd[pix - 1] = truncated(data[ndata == pix])

        avg_rms = self.figure_check_combination(
            rnimage, rmean, rstd, step=step)

        # Fake sky error image
        self.figure_simple_image(ndata, title='Number of frames combined')

        # Create fake error frame
        mask = (ndata <= 0)
        ndata[mask] = 1
        fake = numpy.where(
            mask, 0.0, numpy.random.normal(avg_rms / numpy.sqrt(ndata)))
        ndata[mask] = 0
        self.figure_simple_image(fake, title='Fake sky error image')
        # store fake image
        fits.writeto('fake_sky_rms_i%0d.fits' % step, fake, overwrite=True)

    def figure_check_combination(self, rnimage, rmean, rstd, step=0):
        self._figure.clf()
        self._figure.subplots_adjust(hspace=0.001)

        ax1 = self._figure.add_subplot(3, 1, 1)
        pred = [rstd[-1] * math.sqrt(rnimage[-1] /
                                     float(npix)) for npix in rnimage]
        ax1.plot(rnimage, rstd, 'g*', rnimage, pred, 'y-')
        ax1.set_title("")
        ax1.set_ylabel('measured sky rms')

        ax2 = self._figure.add_subplot(3, 1, 2, sharex=ax1)
        pred = [val * math.sqrt(npix) for val, npix in zip(rstd, rnimage)]
        avg_rms = sum(pred) / len(pred)
        ax2.plot(
            rnimage, pred, 'r*', [rnimage[0], rnimage[-1]], [avg_rms, avg_rms])
        ax2.set_ylabel('scaled sky rms')

        ax3 = self._figure.add_subplot(3, 1, 3, sharex=ax1)
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

    def figure_median_background(self, scales, step=0):
        # FIXME: plotting
        self._figure.clf()
        ax = self._figure.add_subplot(1, 1, 1)
        ax.plot(scales, 'r*')
        ax.set_title("")
        ax.set_xlabel('Image number')
        ax.set_ylabel('Median')
        self._figure.canvas.draw()
        self._figure.savefig('figure-median-sky-background_i%01d.png' % step)

    def figure_image(self, thedata, image):
        ax = self._figure.gca()
        image_axes, = ax.get_images()
        image_axes.set_data(thedata)

        # Create normalizer object
        interval = PercentileInterval(50.)
        z1, z2 = interval.get_limits(thedata)
        norm = ImageNormalize(vmin=z1, vmax=z2, stretch=SqrtStretch())
        image_axes.set_clim(z1, z2)
        image_axes.set_norm(norm)
        clim = image_axes.get_clim()
        ax.set_title('%s, bg=%g fg=%g, linscale' %
                     (image.lastname, clim[0], clim[1]))
        self._figure.canvas.draw()

    def create_mask_single(self, frame, seeing_fwhm, step=0):
        #
        # remove_border = True

        # sextractor takes care of bad pixels
        sex = SExtractor()
        sex.config['CHECKIMAGE_TYPE'] = "SEGMENTATION"
        sex.config["CHECKIMAGE_NAME"] = name_object_mask(frame.baselabel, step)

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

        filename = frame.lastname

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
        self._figure.savefig('figure-sky-segmentation-overlay_%01d.png' % step)

        self.figure_fwhm_histogram(fwhms, step=step)

        # mode with an histogram
        hist, edges = numpy.histogram(fwhms, 50)
        idx = hist.argmax()

        seeing_fwhm = 0.5 * (edges[idx] + edges[idx + 1])
        if seeing_fwhm <= 0:
            _logger.warning(
                'Seeing FHWM %f pixels is negative, reseting', seeing_fwhm)
            seeing_fwhm = None
        else:
            _logger.info('Seeing FHWM %f pixels (%f arcseconds)',
                         seeing_fwhm, seeing_fwhm * sex.config['PIXEL_SCALE'])
        name_segmask(step)
        _logger.info('Step %d, create object mask %s', step,  frame.objmask)
        frame.objmask = name_object_mask(frame.baselabel, step)
        frame.objmask_data = None
        return frame, seeing_fwhm

    def create_mask(self, sf_data, seeing_fwhm, step=0):
        # FIXME more plots
        self.figure_final_before_s(sf_data[0])

        #
        remove_border = True

        # sextractor takes care of bad pixels


        #if seeing_fwhm is not None and seeing_fwhm > 0:
        #    sex.config['SEEING_FWHM'] = seeing_fwhm * sex.config['PIXEL_SCALE']


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
            lower = sf_data[2].max() // 10
            border = (wm < lower)
            fits.writeto(weigthmap, border.astype('uint8'), overwrite=True)

            #sex.config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            # FIXME: this is a magic number
            # sex.config['WEIGHT_THRESH'] = 50
            #sex.config['WEIGHT_IMAGE'] = weigthmap
        else:
            border=None

        filename = 'result_i%0d.fits' % (step)

        # Lauch SExtractor on a FITS file
        #sex.run(filename)

        data_res = fits.getdata(filename)
        data_res = data_res.byteswap().newbyteorder()
        bkg = sep.Background(data_res)
        data_sub = data_res - bkg

        _logger.info('Runing source extraction tor in %s', filename)
        objects, objmask = sep.extract(data_sub, 1.5, err=bkg.globalrms,
                              mask=border, segmentation_map=True)
        fits.writeto(name_segmask(step), objmask, overwrite=True)

        # # Plot objects
        # # FIXME, plot sextractor objects on top of image
        # patches = []
        # fwhms = []
        # nfirst = 0
        # catalog_f = sopen(sex.config['CATALOG_NAME'])
        # try:
        #     star = catalog_f.readline()
        #     while star:
        #         flags = star['FLAGS']
        #         # ignoring those objects with corrupted apertures
        #         if flags & sexcatalog.CORRUPTED_APER:
        #             star = catalog_f.readline()
        #             continue
        #         center = (star['X_IMAGE'], star['Y_IMAGE'])
        #         wd = 10 * star['A_IMAGE']
        #         hd = 10 * star['B_IMAGE']
        #         color = 'red'
        #         e = Ellipse(center, wd, hd, star['THETA_IMAGE'], color=color)
        #         patches.append(e)
        #         fwhms.append(star['FWHM_IMAGE'])
        #         nfirst += 1
        #         # FIXME Plot a ellipse
        #         star = catalog_f.readline()
        # finally:
        #     catalog_f.close()
        #
        # p = PatchCollection(patches, alpha=0.4)
        # ax = self._figure.gca()
        # ax.add_collection(p)
        # self._figure.canvas.draw()
        # self._figure.savefig('figure-segmentation-overlay_%01d.png' % step)
        #
        # self.figure_fwhm_histogram(fwhms, step=step)
        #
        # # mode with an histogram
        # hist, edges = numpy.histogram(fwhms, 50)
        # idx = hist.argmax()
        #
        # seeing_fwhm = 0.5 * (edges[idx] + edges[idx + 1])
        # if seeing_fwhm <= 0:
        #     _logger.warning(
        #         'Seeing FHWM %f pixels is negative, reseting', seeing_fwhm)
        #     seeing_fwhm = None
        # else:
        #     _logger.info('Seeing FHWM %f pixels (%f arcseconds)',
        #                  seeing_fwhm, seeing_fwhm * sex.config['PIXEL_SCALE'])
        # objmask = fits.getdata(name_segmask(step))

        return objmask, seeing_fwhm

    def figure_final_before_s(self, data):
        self._figure.clf()
        ax = self._figure.add_subplot(111)
        cmap = mpl.cm.get_cmap('gray')
        # norm = mpl.colors.LogNorm()
        ax.set_title('Result image')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        interval = PercentileInterval(50.)
        z1, z2 = interval.get_limits(data)
        norm = ImageNormalize(vmin=z1, vmax=z2, stretch=SqrtStretch())
        ax.imshow(data, cmap=cmap, clim=(z1, z2), norm=norm)
        self._figure.canvas.draw()

    def figure_fwhm_histogram(self, fwhms, step=0):
        self._figure.clf()
        ax = self._figure.add_subplot(111)
        ax.set_title('FWHM of objects')
        ax.hist(fwhms, 50, normed=1, facecolor='g', alpha=0.75)
        self._figure.canvas.draw()
        self._figure.savefig('figure-fwhm-histogram_i%01d.png' % step)
