#
# Copyright 2008-2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Recipe for the reduction of imaging mode observations."""

import datetime
import logging
import os
import shutil
import uuid

import numpy
import sep
from astropy.io import fits
from scipy import interpolate
from scipy.spatial import KDTree as KDTree

from numina.core import Result, Requirement, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.core.query import ResultOf
import numina.array as narray
import numina.array.utils as nautils
import numina.array.combine as nacom
import numina.frame.combine as nfcom
from numina.util.context import manage_fits
from numina.frame import resize_fits, custom_region_to_str

import emirdrp.requirements as reqs
import emirdrp.products as prods
from emirdrp.processing.wcs import offsets_from_wcs_imgs
from emirdrp.processing.corr import offsets_from_crosscor_regions
from emirdrp.core.recipe import EmirRecipe

from .naming import (name_redimensioned_frames, name_object_mask,
                     name_skybackground)
from .naming import name_skybackgroundmask, name_skysub_proc, name_skyflat
from .naming import name_skyflat_proc, name_segmask

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2


class ImageInfo(object):
    def __init__(self, origin):
        self.origin = origin
        self.mask = None
        self.metadata = {}
        self.objmask_data = None
        self.valid_target = False
        self.valid_sky = False
        self.label = ""
        self.itype = 'SKY'
        self.valid_region = Ellipsis
        self.resized_base = ""
        self.lastname = ""
        self.flat_corrected = ""

    def __str__(self):
        output = ''
        for item in vars(self):
            value = getattr(self, item)
            output += '\n- {}: {}'.format(item, value)
            if hasattr(value, 'shape'):
                output += ' --> shape {}'.format(value.shape)
        return output


class FullDitheredImagesRecipe(EmirRecipe):
    """Recipe for the reduction of imaging mode observations.

    Recipe to reduce observations obtained in imaging mode, considering
    different possibilities depending on the size of the offsets
    between individual images.
    In particular, the following observing modes are considered: stare imaging,
    nodded beamswitched imaging, and dithered imaging.

    A critical piece of information here is a table that clearly specifies
    which images can be labeled as *science*, and which ones as *sky*.
    Note that some images are used both as *science* and *sky*
    (when the size of the targets is small compared to the offsets).

    **Observing modes:**

     * StareImage
     * Nodded/Beam-switched images
     * Dithered images


    **Inputs:**

     * Science frames + [Sky Frames]
     * Observing mode name: **stare image**, **nodded beamswitched image**,
       or **dithered imaging**
     * A table relating each science image with its sky image(s) (TBD if
       it's in the FITS header and/or in other format)
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

     * Image with three extensions: final image scaled to the individual
       exposure time, variance  and exposure time map OR number of images
       combined (TBD)

    **Procedure:**

    Images are corrected from dark, non-linearity and flat. Then, an iterative
    process starts:

     * Sky is computed from each frame, using the list of sky images of each
       science frame. The objects are avoided using a mask (from the second
       iteration on).

     * The relative offsets are the nominal from the telescope. From the second
       iteration on, we refine them using objects of appropriate brightness
       (not too bright, not to faint).

     * We combine the sky-subtracted images, output is: a new image, a variance
       image and a exposure map/number of images used map.

     * An object mask is generated.

     * We recompute the sky map, using the object mask as an additional input.
       From here we iterate (typically 4 times).

     * Finally, the images are corrected from atmospheric extinction and flux
       calibrated.

     * A preliminary astrometric calibration can always be used (using
       the central coordinates of the pointing and the plate scale
       in the detector).
       A better calibration might be computed using available stars (TBD).

    """

    logger = logging.getLogger(__name__)

    obresult = ObservationResultRequirement(
        query_opts=ResultOf('result_image', node='children')
    )

    master_bpm = reqs.MasterBadPixelMaskRequirement()

    offsets = Requirement(
        prods.CoordinateList2DType,
        'List of pairs of offsets',
        optional=True
    )
    refine_offsets = Parameter(False, 'Refine offsets by cross-correlation')
    iterations = Parameter(0, 'Iterations of the recipe')
    extinction = Parameter(0.0, 'Mean atmospheric extinction')

    method = Parameter(
        'sigmaclip',
        description='Combination method',
        choices=['mean', 'median', 'sigmaclip']
    )
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method',
        optional=True
    )

    sky_images = Parameter(
        0, 'Images used to estimate the '
           'background before and after current image')

    sky_images_sep_time = Parameter(
        10, 'Maximum time interval between target and sky images [minutes]'
    )

    nside_adhoc_sky_correction = Parameter(
        0, 'Ad hoc sky correction (number of subintervals in each quadrant)'
    )

    result_image = Result(prods.ProcessedImage)
    result_sky = Result(prods.ProcessedImage, optional=True)

    def run(self, rinput):

        target_is_sky = True
        obresult = rinput.obresult
        sky_images = rinput.sky_images
        sky_images_sep_time = rinput.sky_images_sep_time
        baseshape = (EMIR_NAXIS2, EMIR_NAXIS1)
        user_offsets = rinput.offsets
        extinction = rinput.extinction
        nside_adhoc_sky_correction = rinput.nside_adhoc_sky_correction

        # protections
        if rinput.iterations == 0 and sky_images != 0:
            raise ValueError('sky_images: {} not compatible with iterations: {}'.format(
                sky_images, rinput.iterations
            ))

        if rinput.iterations > 0 and sky_images == 0:
            raise ValueError('iterations != 0 requires sky_images > 0')

        # check combination method
        if rinput.method != 'sigmaclip':
            if rinput.method_kwargs != {}:
                raise ValueError('Unexpected method_kwargs={}'.format(
                    rinput.method_kwargs))
        # combination method and arguments
        method = getattr(nacom, rinput.method)
        method_kwargs = rinput.method_kwargs

        images_info = self.initial_classification(obresult, target_is_sky)

        # Resizing target frames
        target_info = [iinfo for iinfo in images_info if iinfo.valid_target]
        finalshape, offsetsp, refpix, offset_fc0 = self.compute_size(
            target_info, baseshape, user_offsets
        )

        self.resize(target_info, baseshape, offsetsp, finalshape)

        step = 0

        result = self.process_basic(
            images_info,
            step=step,
            target_is_sky=target_is_sky,
            extinction=extinction,
            method=method,
            method_kwargs=method_kwargs
        )

        if rinput.refine_offsets:
            self.logger.debug("Compute cross-correlation of images")
            # regions_c = self.compute_regions(finalshape, box=200, corners=True)

            # Regions from bright objects
            regions_c = self.compute_regions_from_objs(
                step, result[0].data, finalshape, box=40
            )

            try:

                offsets_xy_c = self.compute_offset_xy_crosscor_regions(
                    images_info, regions_c, refine=True, tol=1
                )
                #
                # Combined offsets
                # Offsets in numpy order, swaping
                offset_xy0 = numpy.fliplr(offset_fc0)
                offsets_xy_t = offset_xy0 - offsets_xy_c
                offsets_fc = numpy.fliplr(offsets_xy_t)
                offsets_fc_t = numpy.round(offsets_fc).astype('int')
                self.logger.debug('Total offsets:\n%s', offsets_xy_t)
                self.logger.info('Computing relative offsets from cross-corr')
                finalshape2, offsetsp2 = narray.combine_shape(
                    baseshape, offsets_fc_t
                )
                #
                self.logger.debug("Relative offsetsp (crosscorr):\n%s",
                                  offsetsp2)
                self.logger.info('Shape of resized array (crosscorr) is '
                                 '(NAXIS2, NAXIS1) = %s', finalshape2)

                # Resizing target imgs
                self.logger.debug("Resize to final offsets")
                self.resize(target_info, baseshape, offsetsp2, finalshape2)
                result = self.process_basic(
                    images_info,
                    step=step,
                    target_is_sky=target_is_sky,
                    extinction=extinction,
                    method=method,
                    method_kwargs=method_kwargs
                )

            except Exception as error:
                self.logger.warning('Error during cross-correlation, %s',
                                    error)

        step = 1

        while step <= rinput.iterations:
            result = self.process_advanced(
                images_info, result, step, target_is_sky,
                maxsep=sky_images_sep_time, nframes=sky_images,
                extinction=extinction,
                method=method, method_kwargs=method_kwargs,
                nside_adhoc_sky_correction=nside_adhoc_sky_correction
            )
            step += 1

        return self.create_result(result_image=result)

    def compute_offset_xy_crosscor_regions(self, iinfo, regions, refine=False,
                                           tol=0.5):

        names = [frame.lastname for frame in iinfo]
        with manage_fits(names) as imgs:
            arrs = [img[0].data for img in imgs]
            offsets_xy = offsets_from_crosscor_regions(
                arrs, regions,
                refine=refine, order='xy', tol=tol
            )
            self.logger.debug("offsets_xy cross-corr:\n%s", offsets_xy)
        return offsets_xy

    def compute_size(self, images_info, baseshape, user_offsets=None):

        # Reference pixel in the center of the frame
        refpix = numpy.array([[baseshape[0] / 2.0, baseshape[1] / 2.0]])

        target_info = [iinfo for iinfo in images_info if iinfo.valid_target]

        if user_offsets is not None:
            self.logger.info('Using offsets from parameters')
            base_ref = numpy.asarray(user_offsets)
            list_of_offsets = -(base_ref - base_ref[0])
        else:
            self.logger.info('Computing offsets from WCS information')
            with manage_fits(img.origin for img in target_info) as images:
                list_of_offsets = offsets_from_wcs_imgs(images, refpix)

        # FIXME: I am using offsets in row/columns
        # the values are provided in XY so flip-lr
        list_of_offsets = numpy.fliplr(list_of_offsets)

        # Insert pixel offsets between frames
        for iinfo, off in zip(target_info, list_of_offsets):
            # Insert pixel offsets between frames
            iinfo.pix_offset = off

            self.logger.debug('Frame %s, offset=%s',
                              iinfo.label, off)

        self.logger.info('Computing relative offsets')
        offsets = [iinfo.pix_offset for iinfo in target_info]
        offsets = numpy.round(offsets).astype('int')

        finalshape, offsetsp = narray.combine_shape(baseshape, offsets)
        self.logger.debug("Relative offsetsp:\n%s", offsetsp)
        self.logger.info('Shape of resized array is (NAXIS2, NAXIS1) = %s',
                         finalshape)
        return finalshape, offsetsp, refpix, list_of_offsets

    def process_basic(self, images_info, step=None,
                      target_is_sky=True,
                      extinction=0.0,
                      method=None, method_kwargs=None):

        target_info = [iinfo for iinfo in images_info if iinfo.valid_target]
        sky_info = [iinfo for iinfo in images_info if iinfo.valid_sky]

        self.logger.info("Step %d, SF: compute superflat", step)
        sf_arr = self.compute_superflat(
            images_info,
            method=method, method_kwargs=method_kwargs
        )

        # Apply superflat
        self.logger.info("Step %d, SF: apply superflat", step)
        for iinfo in images_info:
            self.correct_superflat(iinfo, sf_arr, step=step, save=True)

        self.logger.info('Simple sky correction')
        if target_is_sky:
            # Each frame is the closest sky frame available
            for iinfo in images_info:
                self.compute_simple_sky_for_frame(iinfo, iinfo)
        else:
            # Not implemented
            self.compute_simple_sky(target_info, sky_info)

        # Combining the frames
        self.logger.info("Step %d, Combining target frames", step)
        result = self.combine_frames(target_info, extinction=extinction,
                                     method=method, method_kwargs=method_kwargs)
        self.logger.info('Step %d, finished', step)

        return result

    def process_advanced(self, images_info, result, step, target_is_sky=True,
                         maxsep=5.0, nframes=6, extinction=0,
                         method=None, method_kwargs=None,
                         nside_adhoc_sky_correction=0):

        seeing_fwhm = None
        baseshape = (EMIR_NAXIS2, EMIR_NAXIS1)
        target_info = [iinfo for iinfo in images_info if iinfo.valid_target]
        sky_info = [iinfo for iinfo in images_info if iinfo.valid_sky]
        self.logger.info('Step %d, generating segmentation image', step)

        objmask, seeing_fwhm = self.create_mask(
            result, seeing_fwhm, step=step)

        for frame in target_info:
            frame.objmask = name_object_mask(frame.label, step)
            self.logger.info(
                'Step %d, create object mask %s', step, frame.objmask)
            frame.objmask_data = objmask[frame.valid_region]
            fits.writeto(
                frame.objmask, frame.objmask_data, overwrite=True)

        if not target_is_sky:
            # Empty object mask for sky frames
            bogus_objmask = numpy.zeros(baseshape, dtype='uint8')

            for frame in sky_info:
                frame.objmask_data = bogus_objmask

        self.logger.info("Step %d, SF: compute superflat", step)
        sf_arr = self.compute_superflat(sky_info, segmask=objmask, step=step,
                                        method=method, method_kwargs=method_kwargs)

        # Apply superflat
        self.logger.info("Step %d, SF: apply superflat", step)
        for iinfo in images_info:
            self.correct_superflat(iinfo, sf_arr, step=step, save=True)

        self.logger.info('Step %d, advanced sky correction (SC)', step)
        self.compute_advanced_sky(
            target_info, objmask,
            skyframes=sky_info,
            target_is_sky=target_is_sky,
            maxsep=maxsep,
            nframes=nframes,
            step=step,
            method=method,
            method_kwargs=method_kwargs,
            nside_adhoc_sky_correction=nside_adhoc_sky_correction
        )

        # Combining the images
        self.logger.info("Step %d, Combining the images", step)
        # FIXME: only for science
        result = self.combine_frames(target_info, extinction, step=step,
                                     method=method, method_kwargs=method_kwargs)
        return result

    def compute_simple_sky_for_frame(self, frame, skyframe, step=0, save=True):
        self.logger.info('Correcting sky in frame.....: %s', frame.lastname)
        self.logger.info('with sky computed from frame: %s', skyframe.lastname)

        if hasattr(skyframe, 'median_sky'):
            sky = skyframe.median_sky
        else:

            with fits.open(skyframe.lastname, mode='readonly') as hdulist:
                data = hdulist['primary'].data
                valid = data[frame.valid_region]

                if skyframe.objmask_data is not None:
                    self.logger.debug('object mask defined')
                    msk = frame.objmask_data
                    sky = numpy.median(valid[msk == 0])
                else:
                    self.logger.debug('object mask empty')
                    sky = numpy.median(valid)

            self.logger.debug('median sky value is %f', sky)
            skyframe.median_sky = sky

        dst = name_skysub_proc(frame.label, step)
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
            self.logger.info('Sky-subtrated image in frame: %s',
                             frame.lastname)
            self.logger.info('---')

    def compute_simple_sky(self, frame, skyframe, step=0, save=True):
        raise NotImplementedError

    def correct_superflat(self, frame, fitted, step=0, save=True):

        frame.flat_corrected = name_skyflat_proc(frame.label, step)
        if save:
            shutil.copyfile(frame.resized_base, frame.flat_corrected)
        else:
            os.rename(frame.resized_base, frame.flat_corrected)

        self.logger.info("Step %d, SF: apply superflat to frame %s",
                         step, frame.flat_corrected)
        with fits.open(frame.flat_corrected, mode='update') as hdulist:
            data = hdulist['primary'].data
            datar = data[frame.valid_region]
            data[frame.valid_region] = narray.correct_flatfield(datar, fitted)

            frame.lastname = frame.flat_corrected

    def initial_classification(self, obresult, target_is_sky=False):
        """Classify input frames, """
        # lists of targets and sky frames

        with obresult.frames[0].open() as baseimg:
            # Initial checks
            has_bpm_ext = 'BPM' in baseimg
            self.logger.info('images have BPM extension: %s', has_bpm_ext)

        images_info = []
        for f in obresult.frames:
            with f.open() as img:
                # Getting some metadata from FITS header
                hdr = img[0].header

                iinfo = ImageInfo(f)

                finfo = {}
                iinfo.metadata = finfo

                finfo['uuid'] = hdr['UUID']
                finfo['exposure'] = hdr['EXPTIME']
                # frame.baseshape = get_image_shape(hdr)
                finfo['airmass'] = hdr['airmass']
                finfo['mjd'] = hdr['tstamp']

                iinfo.label = 'result_image_{}'.format(finfo['uuid'])
                iinfo.mask = nfcom.Extension("BPM")
                # Insert pixel offsets between frames
                iinfo.objmask_data = None
                iinfo.valid_target = False
                iinfo.valid_sky = False

                # FIXME: hardcode itype for the moment
                iinfo.itype = 'TARGET'
                if iinfo.itype == 'TARGET':
                    iinfo.valid_target = True
                    # targetframes.append(iinfo)
                    if target_is_sky:
                        iinfo.valid_sky = True
                        # skyframes.append(iinfo)
                if iinfo.itype == 'SKY':
                    iinfo.valid_sky = True
                    # skyframes.append(iinfo)
                images_info.append(iinfo)

        return images_info

    def compute_superflat(self, images_info, segmask=None, step=0,
                          method=None, method_kwargs=None):

        self.logger.info("Step %d, SF: combining the frames without offsets",
                         step)

        base_imgs = [img.resized_base for img in images_info]
        with manage_fits(base_imgs) as imgs:

            data = []
            masks = []

            for img, img_info in zip(imgs, images_info):
                self.logger.debug('Step %d, opening resized frame %s',
                                  step, img_info.resized_base)
                data.append(img['primary'].data[img_info.valid_region])

            scales = [numpy.median(d) for d in data]

            if segmask is not None:
                masks = [segmask[frame.valid_region] for frame in images_info]
            else:
                for frame in images_info:
                    self.logger.debug('Step %d, opening resized mask  %s',
                                      step, frame.resized_mask)
                    hdulist = fits.open(
                         frame.resized_mask, memmap=True, mode='readonly')
                    # filelist.append(hdulist)
                    masks.append(hdulist['primary'].data[frame.valid_region])
                masks = None

            self.logger.debug("Step %d, combining %d frames using '%s'",
                              step, len(data), method.__name__)
            sf_data, _sf_var, sf_num = method(
                data, masks, scales=scales, dtype='float32', **method_kwargs
            )

        # Normalize, flat has mean = 1
        sf_data[sf_data == 0] = 1e-5
        sf_data /= sf_data.mean()
        # sf_data[sf_data <= 0] = 1.0

        # Auxiliary data
        sfhdu = fits.PrimaryHDU(sf_data)
        self.save_intermediate_img(sfhdu, name_skyflat('comb', step))
        return sf_data

    '''
    def compute_sky_advanced(self, data_hdul, omasks, base_header, use_errors):
        method = narray.combine.mean

        self.logger.info("recombine images with segmentation mask using '%s'", method.__name__)
        sky_data = method([m[0].data for m in data_hdul], masks=omasks, dtype='float32')

        hdu = fits.PrimaryHDU(sky_data[0], header=base_header)
        points_no_data = (sky_data[2] == 0).sum()

        self.logger.debug('update created sky image result header')
        skyid = str(uuid.uuid1())
        hdu.header['UUID'] = skyid
        hdu.header['history'] = "Combined {} images using '{}'".format(
            len(data_hdul),
            method.__name__
        )
        hdu.header['history'] = 'Combination time {}'.format(
            datetime.datetime.utcnow().isoformat()
        )
        for img in data_hdul:
            hdu.header['history'] = "Image {}".format(self.datamodel.get_imgid(img))

        msg = "missing pixels, total: {}, fraction: {:3.1f}".format(
            points_no_data,
            points_no_data / sky_data[2].size
        )
        hdu.header['history'] = msg
        self.logger.debug(msg)

        if use_errors:
            varhdu = fits.ImageHDU(sky_data[1], name='VARIANCE')
            num = fits.ImageHDU(sky_data[2], name='MAP')
            sky_result = fits.HDUList([hdu, varhdu, num])
        else:
            sky_result = fits.HDUList([hdu])

        return sky_result
    '''

    def combine_frames(self, frames, extinction, out=None, step=0,
                       method=None, method_kwargs=None):
        self.logger.debug('Step %d, opening sky-subtracted frames', step)

        def fits_open(name):
            """Open FITS with memmap in readonly mode"""
            return fits.open(name, mode='readonly', memmap=True)

        frameslll = [fits_open(frame.lastname)
                     for frame in frames if frame.valid_target]
        self.logger.debug('Step %d, opening mask frames', step)
        mskslll = [fits_open(frame.resized_mask)
                   for frame in frames if frame.valid_target]

        self.logger.debug("Step %d, combining %d frames using '%s'",
                          step, len(frameslll), method.__name__)
        try:
            extinc = [pow(10, -0.4 * frame.metadata['airmass'] * extinction)
                      for frame in frames if frame.valid_target]
            data = [i['primary'].data for i in frameslll]
            masks = [i['primary'].data for i in mskslll]
            headers = [i['primary'].header for i in frameslll]

            out = method(data, masks, scales=extinc, dtype='float32', out=out,
                         **method_kwargs)

            base_header = headers[0]
            hdu = fits.PrimaryHDU(out[0], header=base_header)
            hdu.header['history'] = "Combined %d images using '%s'" % (len(frameslll), method.__name__)
            hdu.header['history'] = 'Combination time {}'.format(datetime.datetime.utcnow().isoformat())
            for img in frameslll:
                hdu.header['history'] = "Image {}".format(img[0].header['uuid'])
            prevnum = base_header.get('NUM-NCOM', 1)
            hdu.header['NUM-NCOM'] = prevnum * len(frameslll)
            hdu.header['NUMRNAM'] = 'FullDitheredImagesRecipe'
            hdu.header['UUID'] = str(uuid.uuid1())
            hdu.header['OBSMODE'] = 'FULL_DITHERED_IMAGE'
            # Headers of last image
            hdu.header['TSUTC2'] = headers[-1]['TSUTC2']

            varhdu = fits.ImageHDU(out[1], name='VARIANCE')
            num = fits.ImageHDU(out[2].astype('uint8'), name='MAP')

            result = fits.HDUList([hdu, varhdu, num])
            # saving the three extensions
            fits.writeto('result_i%0d.fits' % step, out[0], overwrite=True)
            fits.writeto('result_i%0d_var.fits' % step, out[1], overwrite=True)
            fits.writeto('result_i%0d_npix.fits' % step, out[2], overwrite=True)

            result.writeto('result_i%0d_full.fits' % step, overwrite=True)
            return result

        finally:
            self.logger.debug('Step %d, closing sky-subtracted frames', step)
            for f in frameslll:
                f.close()
            self.logger.debug('Step %d, closing mask frames', step)
            for f in mskslll:
                f.close()

    def resize(self, frames, shape, offsetsp, finalshape, window=None,
               scale=1, step=0):
        self.logger.info('Resizing frames and masks')
        self.logger.debug('shape, finalshape (NAXIS2, NAXIS1) = %s --> %s',
                          shape, finalshape)
        for frame, rel_offset in zip(frames, offsetsp):
            if frame.valid_target:
                region, _ = narray.subarray_match(finalshape, rel_offset, shape)
                # Valid region
                frame.valid_region = region
                # Relative offset
                frame.rel_offset = rel_offset
                # names of frame and mask
                framen, maskn = name_redimensioned_frames(
                    frame.label, step)
                frame.resized_base = framen
                frame.resized_mask = maskn
                self.logger.debug('%s', frame.label)
                self.logger.debug('valid region is %s, relative offset is %s',
                                  custom_region_to_str(region), rel_offset)
                self.resize_frame_and_mask(
                    frame, finalshape, framen, maskn, window, scale)
                self.logger.debug('---')

    def resize_frame_and_mask(self, frame, finalshape,
                              framen, maskn, window, scale):
        self.logger.info('Resizing frame %s', frame.label)
        with frame.origin.open() as hdul:
            baseshape = hdul[0].data.shape

            # FIXME: Resize_fits saves the resized image in framen
            resize_fits(hdul, framen, finalshape, frame.valid_region,
                        window=window, scale=scale, dtype='float32')

        self.logger.info('Resizing mask  %s', frame.label)
        # We don't conserve the sum of the values of the frame here, just
        # expand the mask

        if frame.mask is None:
            self.logger.warning('BPM missing, use zeros instead')
            false_mask = numpy.zeros(baseshape, dtype='int16')
            hdum = fits.HDUList(fits.PrimaryHDU(false_mask))
            frame.mask = hdum  # DataFrame(frame=hdum)
        elif isinstance(frame.mask, nfcom.Extension):
            ename = frame.mask.name
            with frame.origin.open() as hdul:
                frame.mask = fits.HDUList(hdul[ename].copy())

        resize_fits(frame.mask, maskn, finalshape, frame.valid_region,
                    fill=1, window=window, scale=scale, conserve=False)

    def create_mask(self, img, seeing_fwhm, step=0):

        #
        remove_border = True

        # sextractor takes care of bad pixels

        # if seeing_fwhm is not None and seeing_fwhm > 0:
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
            wm = img[2].data.copy()
            # Dont search objects where nimages < lower
            # FIXME: this is a magic number
            # We ignore objects in regions where we have less
            # than 10% of the images
            lower = wm.max() // 10
            border = (wm < lower)
            fits.writeto(weigthmap, border.astype('uint8'), overwrite=True)

            # sex.config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            # FIXME: this is a magic number
            # sex.config['WEIGHT_THRESH'] = 50
            # sex.config['WEIGHT_IMAGE'] = weigthmap
        else:
            border = None

        data_res = img[0].data
        bkg = sep.Background(data_res)
        data_sub = data_res - bkg

        self.logger.info('Runing source extraction in previous result')
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

    def compute_advanced_sky(self, targetframes, objmask,
                             skyframes=None, target_is_sky=False,
                             maxsep=5.0,
                             nframes=10,
                             step=0,
                             method=None, method_kwargs=None,
                             nside_adhoc_sky_correction=0):

        if target_is_sky:
            skyframes = targetframes
            # Each frame is its closest sky frame
            nframes += 1
        elif skyframes is None:
            raise ValueError('skyframes not defined')

        # build kdtree
        sarray = numpy.array([frame.metadata['mjd'] for frame in skyframes])
        # shape must be (n, 1)
        sarray = numpy.expand_dims(sarray, axis=1)

        # query
        tarray = numpy.array([frame.metadata['mjd'] for frame in targetframes])
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
                self.logger.info("Step %d, SC: computing advanced sky for %s using '%s'",
                                 step, tf.label, method.__name__)
                # filter(lambda x: x < nsky, idss)
                locskyframes = []
                for si in idss:
                    if tid == si:
                        # this sky frame it is the current frame, reject
                        continue
                    if si < nsky:
                        self.logger.debug('Step %d, SC: %s is a sky frame',
                                          step, skyframes[si].label)
                        locskyframes.append(skyframes[si])
                self.compute_advanced_sky_for_frame(
                    tf, locskyframes, step=step,
                    method=method, method_kwargs=method_kwargs,
                    nside_adhoc_sky_correction=nside_adhoc_sky_correction
                )
            except IndexError:
                self.logger.error('No sky image available for frame %s', tf.lastname)
                raise

    def compute_advanced_sky_for_frame(self, frame, skyframes, step=0,
                                       method=None, method_kwargs=None,
                                       nside_adhoc_sky_correction=0):
        self.logger.info('Correcting sky in frame %s', frame.lastname)
        self.logger.info('with sky computed from frames')
        for i in skyframes:
            self.logger.info('%s', i.flat_corrected)

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
                    self.logger.debug('object mask is shared')
                elif i.objmask is not None:
                    hdulistmask = fits.open(
                        i.objmask, mode='readonly', memmap=True)
                    masks.append(hdulistmask['primary'].data)
                    desc.append(hdulistmask)
                    self.logger.debug('object mask is particular')
                else:
                    self.logger.warn('no object mask for %s', filename)

            self.logger.debug("computing background with %d frames using '%s'",
                              len(data), method.__name__)
            sky, _, num = method(data, masks, scales=scales, **method_kwargs)

            with fits.open(frame.lastname) as hdulist:
                data = hdulist['primary'].data
                valid = data[frame.valid_region]

                if frame.objmask_data is not None:
                    self.logger.debug('object mask defined')
                    msk = frame.objmask_data
                    skymedian = numpy.median(valid[msk == 0])
                else:
                    self.logger.debug('object mask empty')
                    skymedian = numpy.median(valid)
                self.logger.debug('scaling with skymedian %s', skymedian)
            sky *= skymedian

        finally:
            # Closing all FITS files
            for hdl in desc:
                hdl.close()

        if numpy.any(num == 0):
            # We have pixels without
            # sky background information
            self.logger.warning('pixels without sky information when correcting %s',
                                frame.flat_corrected)
            binmask = num == 0
            # FIXME: during development, this is faster
            # sky[binmask] = sky[num != 0].mean()

            # To continue we interpolate over the patches
            narray.fixpix2(sky, binmask, out=sky, iterations=1)

            name = name_skybackgroundmask(frame.label, step)
            fits.writeto(name, binmask.astype('int16'), overwrite=True)

        name_sky = name_skybackground(frame.label, step)
        fits.writeto(name_sky, sky, overwrite=True)

        dst = name_skysub_proc(frame.label, step)
        prev = frame.lastname
        shutil.copyfile(prev, dst)
        frame.lastname = dst

        with fits.open(frame.lastname, mode='update') as hdulist:
            data = hdulist['primary'].data
            valid = data[frame.valid_region]
            valid -= sky
            if nside_adhoc_sky_correction > 0:
                skycorr = self.adhoc_sky_correction(
                    arr=valid,
                    objmask=frame.objmask_data,
                    nside=nside_adhoc_sky_correction
                )
                valid -= skycorr

    def compute_regions_from_objs(self, step, arr, finalshape, box=50,
                                  corners=True):
        regions = []
        # create catalog of objects skipping a border around the image
        catalog, mask = self.create_object_catalog(arr, border=300)

        self.save_intermediate_array(mask, 'objmask_i{}.fits'.format(step))
        # with the catalog, compute the brightest NKEEP objects

        LIMIT_AREA = 5000
        NKEEP = 3
        idx_small = catalog['npix'] < LIMIT_AREA
        objects_small = catalog[idx_small]
        idx_flux = objects_small['flux'].argsort()
        objects_nth = objects_small[idx_flux][-NKEEP:]
        for obj in objects_nth:
            self.logger.debug('ref is (x,y) = (%s, %s)', obj['x'], obj['y'])
            region = nautils.image_box2d(obj['x'], obj['y'], finalshape, (box, box))
            regions.append(region)
        return regions

    def create_object_catalog(self, arr, threshold=3.0, border=0):

        if border > 0:
            wmap = numpy.ones_like(arr)
            wmap[border:-border, border:-border] = 0
        else:
            wmap = None

        bkg = sep.Background(arr)
        data_sub = arr - bkg
        objects, objmask = sep.extract(
            data_sub,
            threshold,
            err=bkg.globalrms * numpy.ones_like(data_sub),
            mask=wmap,
            segmentation_map=True
        )
        return objects, objmask

    def adhoc_sky_correction(self, arr, objmask, nside=10):
        # nside: number of subdivisions in each quadrant

        self.logger.info('computing ad hoc sky correction')

        skyfit = numpy.zeros_like(arr)

        # fit each quadrant separately
        lim = [0, 1024, 2048]
        for i in range(2):
            i1 = lim[i]
            i2 = lim[i + 1]
            for j in range(2):
                j1 = lim[j]
                j2 = lim[j + 1]
                xyfit = []
                zfit = []
                limi = numpy.linspace(i1, i2, nside + 1, dtype=int)
                limj = numpy.linspace(j1, j2, nside + 1, dtype=int)
                for ii in range(nside):
                    ii1 = limi[ii]
                    ii2 = limi[ii+1]
                    for jj in range(nside):
                        jj1 = limj[jj]
                        jj2 = limj[jj+1]
                        tmprect = arr[ii1:ii2, jj1:jj2]
                        tmpmask = objmask[ii1:ii2, jj1:jj2]
                        usefulpix = tmprect[tmpmask == 0]
                        if usefulpix.size > 0:
                            x0 = (jj1 + jj2) / 2.0
                            y0 = (ii1 + ii2) / 2.0
                            z0 = numpy.median(usefulpix)
                            xyfit.append([x0, y0])
                            zfit.append(z0)
                xyfit = numpy.array(xyfit)
                zfit = numpy.array(zfit)
                xgrid, ygrid = numpy.meshgrid(
                    numpy.arange(j1, j2, dtype=float),
                    numpy.arange(i1, i2, dtype=float))
                surface_nearest = interpolate.griddata(
                    xyfit, zfit, (xgrid, ygrid), method='nearest',
                    rescale=True
                )
                surface_cubic = interpolate.griddata(
                    xyfit, zfit, (xgrid, ygrid), method='cubic',
                    fill_value=-1.0E30,
                    rescale=True
                )
                skyfit[i1:i2, j1:j2] = numpy.where(
                    surface_cubic < -1.0E29,
                    surface_nearest,
                    surface_cubic
                )

        # hdu = fits.PrimaryHDU(arr.astype('float32'))
        # hdul = fits.HDUList([hdu])
        # hdul.writeto('xxx1.fits', overwrite=True)
        # hdu = fits.PrimaryHDU(skyfit.astype('float32'))
        # hdul = fits.HDUList([hdu])
        # hdul.writeto('xxx2.fits', overwrite=True)

        return skyfit
