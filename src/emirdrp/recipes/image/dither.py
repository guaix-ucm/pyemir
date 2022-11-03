#
# Copyright 2008-2022 Universidad Complutense de Madrid
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
from scipy.ndimage import median_filter
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
        self.pix_offset = None
        self.rel_offset = None
        self.resized_base = ""
        self.resized_mask = ""
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

    **Note (June 2022)**: the following documentation is not updated. It
    contains the initial goals, but some of them are still not implemented.

    --- Not updated ---

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
        query_opts=ResultOf('reduced_image', node='children')
    )

    master_bpm = reqs.MasterBadPixelMaskRequirement()

    offsets = Requirement(
        prods.CoordinateList2DType,
        'List of pairs of offsets',
        optional=True
    )
    refine_offsets = Parameter(False, 'Refine offsets by cross-correlation')
    iterations = Parameter(0, 'Iterations of the recipe')
    fit_doughnut = Parameter(False, 'Fit doughnut-like shape in superflat')
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

    reduced_image = Result(prods.ProcessedImage)
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
        finalshape, offsetsp, offset_fc0 = self.compute_size(
            target_info, baseshape, user_offsets
        )

        self.resize_all(target_info, baseshape, offsetsp, finalshape)

        step = 0

        result = self.process_basic(
            images_info,
            step=step,
            target_is_sky=target_is_sky,
            extinction=extinction,
            method=method,
            method_kwargs=method_kwargs,
            fit_doughnut=rinput.fit_doughnut
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
                self.logger.debug(f"Relative offsetsp (crosscorr):\n{offsetsp2}"),
                self.logger.info(f'Shape of resized array (crosscorr) is (NAXIS2, NAXIS1) = {finalshape2}')

                # Resizing target imgs
                self.logger.debug("Resize to final offsets")
                self.resize_all(target_info, baseshape, offsetsp2, finalshape2)
                result = self.process_basic(
                    images_info,
                    step=step,
                    target_is_sky=target_is_sky,
                    extinction=extinction,
                    method=method,
                    method_kwargs=method_kwargs,
                    fit_doughnut=rinput.fit_doughnut
                )

            except Exception as error:
                self.logger.warning(f'Error during cross-correlation, {error}')

        step = 1

        while step <= rinput.iterations:
            result = self.process_advanced(
                images_info, result, step, target_is_sky,
                maxsep_time=sky_images_sep_time, nframes=sky_images,
                extinction=extinction,
                method=method, method_kwargs=method_kwargs,
                nside_adhoc_sky_correction=nside_adhoc_sky_correction,
                fit_doughnut=rinput.fit_doughnut
            )
            step += 1

        return self.create_result(reduced_image=result)

    def compute_offset_xy_crosscor_regions(self, iinfo, regions, refine=False, tol=0.5):

        names = [frame.lastname for frame in iinfo]
        with manage_fits(names) as imgs:
            arrs = [img[0].data for img in imgs]
            offsets_xy = offsets_from_crosscor_regions(
                arrs, regions,
                refine=refine, order='xy', tol=tol
            )
            self.logger.debug(f"offsets_xy cross-corr:\n{offsets_xy}")
        return offsets_xy

    def compute_size(self, target_info, baseshape, user_offsets=None):
        """Compute relative offsets between images. """

        if user_offsets is not None:
            self.logger.info('Using offsets from parameters')
            base_ref = numpy.asarray(user_offsets)
            list_of_offsets = -(base_ref - base_ref[0])
        else:
            self.logger.info('Computing offsets from WCS information')
            # Reference pixel in the center of the frame
            refpix = numpy.array([[baseshape[0] / 2.0, baseshape[1] / 2.0]])
            with manage_fits(iinfo.origin for iinfo in target_info) as images:
                list_of_offsets = offsets_from_wcs_imgs(images, refpix)

        # the values are provided in XY so flip-lr
        list_of_offsets = numpy.fliplr(list_of_offsets)

        # store pixel offsets between frames
        if len(target_info) != len(list_of_offsets):
            raise ValueError(
                f'The number of offsets ({len(list_of_offsets)})'
                f' does not match the number of images ({len(target_info)})'
            )
        for iinfo, yx_offset in zip(target_info, list_of_offsets):
            # Insert pixel offsets between frames
            iinfo.pix_offset = yx_offset
            self.logger.debug(f'Frame {iinfo.label}, offset [Y, X] = {yx_offset}')

        self.logger.info('Computing relative offsets')
        offsets = [iinfo.pix_offset for iinfo in target_info]
        offsets = numpy.round(offsets).astype('int')

        finalshape, offsetsp = narray.combine_shape(baseshape, offsets)
        self.logger.debug(f"Relative offsetsp [Y, X] (from lower left corner):\n{offsetsp}")
        self.logger.info(f'Shape of resized array is (NAXIS2, NAXIS1) = {finalshape}')
        return finalshape, offsetsp, list_of_offsets

    def process_basic(self, images_info, step=None,
                      target_is_sky=True,
                      extinction=0.0,
                      method=None, method_kwargs=None,
                      fit_doughnut=False):

        target_info = [iinfo for iinfo in images_info if iinfo.valid_target]
        sky_info = [iinfo for iinfo in images_info if iinfo.valid_sky]

        self.logger.info(f"Step {step}, SF: compute superflat")
        sf_arr, doughnut_arr = self.compute_superflat(
            images_info,
            method=method, method_kwargs=method_kwargs,
            fit_doughnut=fit_doughnut
        )

        # Apply superflat
        self.logger.info(f"Step {step}, SF: apply superflat")
        for iinfo in images_info:
            self.correct_superflat(iinfo, sf_arr, step=step, save=True)

        self.logger.info('Simple sky correction')
        self.logger.info('---')
        if target_is_sky:
            # Each frame is the closest sky frame available
            for ii, iinfo in enumerate(images_info):
                self.logger.info(f'image {ii+1} / {len(images_info)}')
                self.compute_simple_sky_for_frame(iinfo, iinfo, doughnut_arr, step)
                self.logger.info('---')
        else:
            # Not implemented
            self.compute_simple_sky(target_info, sky_info, step)

        # Combining the frames
        self.logger.info(f"Step {step}, Combining target frames")
        result = self.combine_frames(target_info, extinction=extinction,
                                     method=method, method_kwargs=method_kwargs)
        self.logger.info(f"Step {step}, finished")
        self.logger.info('---')

        return result

    def process_advanced(self, images_info, result, step, target_is_sky=True,
                         maxsep_time=5.0, nframes=6, extinction=0,
                         method=None, method_kwargs=None,
                         nside_adhoc_sky_correction=0,
                         fit_doughnut=False):

        seeing_fwhm = None
        baseshape = (EMIR_NAXIS2, EMIR_NAXIS1)
        target_info = [iinfo for iinfo in images_info if iinfo.valid_target]
        sky_info = [iinfo for iinfo in images_info if iinfo.valid_sky]
        self.logger.info(f'Step {step}, generating segmentation image')

        # create object mask from combined result
        objmask, seeing_fwhm = self.create_objmask(
            result, seeing_fwhm, step=step)

        for frame in target_info:
            frame.objmask = name_object_mask(frame.label, step)
            self.logger.info(f'Step {step}, create object mask (+ footprint) {frame.objmask}')
            footprint = frame.mask[0].data
            if footprint.shape != baseshape:
                raise ValueError(f'Unexpected footprint.shape: {footprint.shape} != {baseshape}')
            # important: include footprint in objmask
            frame.objmask_data = objmask[frame.valid_region] + footprint
            fits.writeto(
                frame.objmask, frame.objmask_data, overwrite=True)

        if not target_is_sky:
            # Empty object mask for sky frames
            bogus_objmask = numpy.zeros(baseshape, dtype='uint8')

            for frame in sky_info:
                footprint = frame.mask[0].data
                if footprint.shape != baseshape:
                    raise ValueError(f'Unexpected footprint.shape: {footprint.shape} != {baseshape}')
                frame.objmask_data = bogus_objmask + footprint

        self.logger.info(f"Step {step}, SF: compute superflat")
        sf_arr, doughnut_arr = self.compute_superflat(
            sky_info,
            segmask=objmask,
            step=step,
            method=method,
            method_kwargs=method_kwargs,
            fit_doughnut=fit_doughnut
        )

        # Apply superflat
        self.logger.info(f"Step {step}, SF: apply superflat")
        for iinfo in images_info:
            self.correct_superflat(iinfo, sf_arr, step=step, save=True)

        self.logger.info(f'Step {step}, Advanced sky correction (SC)')
        self.compute_advanced_sky(
            target_info, objmask,
            skyframes=sky_info,
            target_is_sky=target_is_sky,
            maxsep_time=maxsep_time,
            nframes=nframes,
            step=step,
            method=method,
            method_kwargs=method_kwargs,
            nside_adhoc_sky_correction=nside_adhoc_sky_correction
        )

        # Combining the images
        self.logger.info('---')
        self.logger.info(f"Step {step}, Combining the images")
        # FIXME: only for science
        result = self.combine_frames(target_info, extinction, step=step,
                                     method=method, method_kwargs=method_kwargs)
        return result

    def compute_simple_sky_for_frame(self, frame, skyframe, doughnut_arr=None, step=0, save=True):
        self.logger.info(f'Correcting sky in frame.....: {frame.lastname}')
        self.logger.info(f'with sky computed from frame: {skyframe.lastname}')

        if hasattr(skyframe, 'median_sky'):
            sky = skyframe.median_sky
            self.logger.debug(f'using previously computed median sky {sky}')
        else:

            with fits.open(skyframe.lastname, mode='readonly') as hdulist:
                data = hdulist['primary'].data
                valid = data[frame.valid_region]

                if skyframe.objmask_data is not None:
                    self.logger.debug('object mask defined (it must include the footprint)')
                    msk = frame.objmask_data
                else:
                    self.logger.debug('object mask empty (using only footprint)')
                    footprint = skyframe.mask[0].data
                    msk = footprint
                sky = numpy.median(valid[msk == 0])

            self.logger.debug(f'median sky value is {sky}')
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
            if doughnut_arr is None:
                valid -= sky
            else:
                # since the doughnut_arr was computed when deriving the superflat,
                # its averaged signal is around 1.0; for that reason here we
                # scale it by the median sky value
                self.logger.debug(f'subtracting (median sky value) * doughnut_fit')
                valid -= sky * doughnut_arr
            self.logger.info(f'Sky-subtrated image in frame: {frame.lastname}')

    def compute_simple_sky(self, frame, skyframe, step=0, save=True):
        raise NotImplementedError

    def correct_superflat(self, frame, fitted, step=0, save=True):

        frame.flat_corrected = name_skyflat_proc(frame.label, step)
        if save:
            shutil.copyfile(frame.resized_base, frame.flat_corrected)
        else:
            os.rename(frame.resized_base, frame.flat_corrected)

        self.logger.info(f"Step {step}, SF: apply superflat, generating {frame.flat_corrected}")
        with fits.open(frame.flat_corrected, mode='update') as hdulist:
            data = hdulist['primary'].data
            datar = data[frame.valid_region]
            # although the superflat contains very small values (1e-5) outside
            # the useful footprint region (in order to avoid division by zero),
            # those pixels in the array datar are zero and the flatfield corrected
            # result is correct
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

                iinfo.label = 'reduced_image_{}'.format(finfo['uuid'])
                iinfo.mask = nfcom.Extension("BPM")
                # Insert pixel offsets between frames
                iinfo.objmask_data = None
                iinfo.valid_target = False
                iinfo.valid_sky = False

                # ToDo: revise this!
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
                          method=None, method_kwargs=None, fit_doughnut=False):

        self.logger.info(f"Step {step}, SF: combining the frames without offsets")

        base_imgs = [img.resized_base for img in images_info]
        with manage_fits(base_imgs) as imgs:

            data = []
            scales = []
            for img, img_info in zip(imgs, images_info):
                # read resized data and select valid_region
                self.logger.debug(f'Step {step}, opening resized frame {img_info.resized_base}')
                tmp_data = img['primary'].data[img_info.valid_region]
                data.append(tmp_data)
                # to compute the proper scale, it is important to skip the
                # masked pixels (and those outside the image footprint);
                # read resized mask and select valid_region
                self.logger.debug(f'Step {step}, opening resized mask  {img_info.resized_mask}')
                with fits.open(img_info.resized_mask, mode='readonly') as hdul:
                    tmp_mask = hdul['primary'].data[img_info.valid_region]
                scales.append(numpy.median(tmp_data[tmp_mask == 0]))

            self.logger.debug(f'Step {step}, scales: {scales}')

            if segmask is not None:
                # segmask contains the object mask (when iterating, step > 0) for
                # the full combined image (final size); for that reason it is important
                # to extract the valid_region for each individual exposure
                masks = [segmask[frame.valid_region] for frame in images_info]
            else:
                '''
                masks = []
                for frame in images_info:
                    self.logger.debug('Step %d, opening resized mask  %s',
                                      step, frame.resized_mask)
                    hdulist = fits.open(
                         frame.resized_mask, memmap=True, mode='readonly')
                    masks.append(hdulist['primary'].data[frame.valid_region])
                '''
                # finally we are not using masks here because the masked
                # pixels in the bad-pixel mask were interpolated with the
                # StareImageRecipe2, whereas the masked pixels corresponding
                # to those outside the footprint of the reprojected images
                # contain zeros (there is no problem averaging them)
                masks = None

            # note that the only masks relevant here are the object masks;
            # the footprint mask is not required because the data outside
            # the image footprint is zero in all the individual exposures to
            # be combined (the computed superflat is initially zero in those
            # pixels outside the footprint region: for that reason these pixels
            # are set to a small, but not zero, value below)
            self.logger.debug(f"Step {step}, combining {len(data)} frames using '{method.__name__}'")
            sf_data, _sf_var, sf_num = method(
                data, masks, scales=scales, dtype='float32', **method_kwargs
            )
            # avoid pixels without flatfield information
            if numpy.any(sf_num == 0):
                self.logger.warning('pixels without flatfield information found: potential problem!')
                self.logger.warning('interpolating missing flatfield pixels using neighbouring data')
                binmask = sf_num == 0
                narray.fixpix2(sf_data, binmask, out=sf_data, iterations=1)

        # Normalize, flat has mean = 1
        mean_signal = numpy.median(sf_data[sf_data > 0])   # avoid region outside footprint
        sf_data /= mean_signal
        # avoid division by zero when using the superflat
        sf_data[sf_data <= 0] = 1E-5
        # sf_data[sf_data <= 0] = 1.0   # do not use this

        # Save superflat
        sfhdu = fits.PrimaryHDU(sf_data)
        tmp_filename = name_skyflat('comb', step)
        self.logger.debug(f"Step {step}, saving {tmp_filename}")
        self.save_intermediate_img(sfhdu, tmp_filename)

        # Compute doughnut fit
        if fit_doughnut:
            self.logger.debug(f"Step {step}, fitting doughnut-like surface")
            doughnut_arr = self.fit_sf_doughnut(image=sf_data, method='linear', fill_with_nearest=True)
            # save doughnut fit
            sfdhdu = fits.PrimaryHDU(doughnut_arr)
            tmp_filename = name_skyflat('doughnut', step)
            self.logger.debug(f"Step {step}, saving {tmp_filename}")
            self.save_intermediate_img(sfdhdu, tmp_filename)
            # recompute superflat
            sf_data /= doughnut_arr
            # save dougnut-corrected superflat
            sfhdu = fits.PrimaryHDU(sf_data)
            tmp_filename = name_skyflat('comb_dc', step)
            self.logger.debug(f"Step {step}, saving {tmp_filename}")
            self.save_intermediate_img(sfhdu, tmp_filename)
        else:
            doughnut_arr = None

        # Return superflat
        return sf_data, doughnut_arr

    def fit_sf_doughnut(self, image, method=None, fill_with_nearest=True):

        if method not in ['test', 'rg_linear', 'nearest', 'linear', 'cubic']:
            raise ValueError(f'Unexpected reconstruction_method={method}')

        # avoid undefined regions (due to image reprojection, for example)
        # note: we have employed 1E-5 in compute_superflat() to avoid
        #       division by zero; for that reason here we use 2E-5
        footprint_useful = (image > 2E-5)

        # image dimensions and center
        naxis2, naxis1 = image.shape
        xc = naxis1 / 2 + 0.5
        yc = naxis2 / 2 + 0.5

        pixel_x = numpy.arange(1, naxis1 + 1)
        pixel_y = numpy.arange(1, naxis2 + 1)

        ix_array = numpy.meshgrid(pixel_x, pixel_y)[0].flatten()
        iy_array = numpy.meshgrid(pixel_x, pixel_y)[1].flatten()

        # polar coordinates
        r = numpy.sqrt((ix_array - xc) ** 2 + (iy_array - yc) ** 2)
        theta = numpy.arctan2(iy_array - yc, ix_array - xc)

        # 2D histogram using r, theta
        nbins_r = 140
        r_bins = numpy.linspace(0, max(r), nbins_r + 1)

        nbins_theta_per_quarter = 90
        nbins_theta = nbins_theta_per_quarter * 4
        theta_bins = numpy.linspace(-180, 180, nbins_theta + 1) * numpy.pi / 180

        hist2d = numpy.zeros((nbins_theta, nbins_r))

        imax = -numpy.ones(nbins_theta, dtype=int)
        for i in range(nbins_r):
            iok = (r >= r_bins[i]) & (r < r_bins[i + 1])
            if numpy.sum(iok) > 0:
                xdum = r[iok]
                ydum = theta[iok]
                zdum = image.flatten()[iok]
                fdum = footprint_useful.flatten()[iok]
                for j in range(nbins_theta):
                    jok = (ydum >= theta_bins[j]) & (ydum < theta_bins[j + 1]) & fdum
                    if numpy.sum(jok) > 0:
                        hist2d[j, i] = numpy.median(zdum[jok])
                        if i > imax[j]:
                            imax[j] = i
                    else:
                        hist2d[j, i] = -1
            else:
                for j in range(nbins_theta):
                    hist2d[j, i] = -1

        # fill rows
        for j in range(nbins_theta):
            for i in range(imax[j] + 1, nbins_r):
                hist2d[j, i] = hist2d[j, imax[j]]

        # for the first radii (with empty data in some bins), average all
        # undefined values with the median within each quadrant
        j1, j2, j3, j4, j5 = numpy.arange(5 * nbins_theta_per_quarter, step=nbins_theta_per_quarter)
        for i in range(nbins_r):
            if len(numpy.argwhere(hist2d[:, i] == -1)) > 0:
                for jj1, jj2 in zip([j1, j2, j3, j4],
                                    [j2, j3, j4, j5]):
                    minicolumn = hist2d[jj1:jj2, i]
                    hist2d[jj1:jj2, i] = numpy.median(minicolumn[minicolumn > -1])

        # smooth result in theta
        hist2d_smooth = median_filter(hist2d, size=(21, 1), mode='wrap')

        # reconstruct doughnut
        if method == 'test':
            # assign to each pixel the value corresponding to hist2d_smooth
            ii = numpy.searchsorted(r_bins, r) - 1
            jj = numpy.searchsorted(theta_bins, theta) - 1
            image2d = hist2d_smooth[jj, ii].reshape(naxis2, naxis1)
        elif method == 'rg_linear':
            # centers of the bins used to compute hist2d_smooth
            r_values = (r_bins[:-1] + r_bins[1:]) / 2
            theta_values = (theta_bins[:-1] + theta_bins[1:]) / 2
            interp = interpolate.RegularGridInterpolator(
                (r_values, theta_values),
                hist2d_smooth.T,
                bounds_error=False
            )
            image2d = interp(numpy.column_stack([r, theta]), method='linear').reshape(naxis2, naxis1)
        else:
            # use scipy.interpolate.gridddata
            zfit = hist2d_smooth.flatten()
            # centers of the bins used to compute hist2d_smooth
            r_values = (r_bins[:-1] + r_bins[1:]) / 2
            theta_values = (theta_bins[:-1] + theta_bins[1:]) / 2
            # cartesian coordinates
            tmp_r = numpy.meshgrid(r_values, theta_values)[0].flatten()
            tmp_theta = numpy.meshgrid(r_values, theta_values)[1].flatten()
            xx = tmp_r * numpy.cos(tmp_theta) + xc
            yy = tmp_r * numpy.sin(tmp_theta) + yc
            # perform fit
            grid_x, grid_y = numpy.mgrid[0:naxis1, 0:naxis2]
            xyfit = numpy.column_stack((xx - 1, yy - 1))
            image2d = interpolate.griddata(xyfit, zfit, (grid_x, grid_y), method=method)
            if fill_with_nearest:
                image2d_fill = interpolate.griddata(xyfit, zfit, (grid_x, grid_y), method='nearest')
                invalid = numpy.isnan(image2d)
                image2d[invalid] = image2d_fill[invalid]
            # important: note that here we are using np.mgrid() instead of np.meshgrid()
            # (see discussion in http://louistiao.me/posts/numpy-mgrid-vs-meshgrid/)
            # and we need to take the transpose
            image2d = image2d.T

        return image2d

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

        # define auxiliary function
        def fits_open(name):
            """Open FITS with memmap in readonly mode"""
            return fits.open(name, mode='readonly', memmap=True)

        self.logger.debug(f'Step {step}, opening sky-subtracted frames')
        frameslll = [fits_open(frame.lastname)
                     for frame in frames if frame.valid_target]

        self.logger.debug(f'Step {step}, opening mask frames')
        mskslll = [fits_open(frame.resized_mask)
                   for frame in frames if frame.valid_target]

        self.logger.debug(f"Step {step}, combining {len(frameslll)} frames using '{method.__name__}'")
        try:
            extinc = [pow(10, -0.4 * frame.metadata['airmass'] * extinction)
                      for frame in frames if frame.valid_target]
            data = [i['primary'].data for i in frameslll]
            masks = [i['primary'].data for i in mskslll]
            headers = [i['primary'].header for i in frameslll]

            # compute combination
            out = method(data, masks, scales=extinc, dtype='float32', out=out,
                         **method_kwargs)

            # update header
            base_header = headers[0]
            hdu = fits.PrimaryHDU(out[0], header=base_header)
            # remove previous HISTORY entries
            # (corresponding to the reduction of the first indidivual exposure)
            self.logger.debug(f"Step {step}, preserving primary header from first exposure")
            self.logger.debug(f"Step {step}, removing HISTORY entries in previous header")
            while 'HISTORY' in hdu.header:
                hdu.header.remove('history')
            # define new HISTORY entries with the id of the combined images
            self.logger.debug(f"Step {step}, updating HISTORY entries with list of individual exposures")
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
            # saving the combined image (result)
            result.writeto('result_i%0d_full.fits' % step, overwrite=True)
            return result

        finally:
            self.logger.debug(f'Step {step}, closing sky-subtracted frames')
            for f in frameslll:
                f.close()
            self.logger.debug(f'Step {step}, closing mask frames')
            for f in mskslll:
                f.close()

    def resize_all(self, target_info, shape, offsetsp, finalshape, window=None,
               scale=1, step=0):
        """Insert each exposure within an array with the final shape """

        self.logger.info('Resizing frames and masks')
        self.logger.debug(f'shape, finalshape (NAXIS2, NAXIS1) = {shape} --> {finalshape}')
        self.logger.debug('---')
        for iframe, (iinfo, rel_offset) in enumerate(zip(target_info, offsetsp)):
            if iinfo.valid_target:
                self.logger.debug(f'image {iframe+1} / {len(offsetsp)}')
                region, _ = narray.subarray_match(finalshape, rel_offset, shape)
                # Valid region
                iinfo.valid_region = region
                # Relative offset
                iinfo.rel_offset = rel_offset
                # names of frame and mask
                frame_name, mask_name = name_redimensioned_frames(iinfo.label, step)
                iinfo.resized_base = frame_name
                iinfo.resized_mask = mask_name
                self.logger.debug(f'{iinfo.label}')
                self.logger.debug(f'valid region (in Python format) is {custom_region_to_str(region)},'
                                  f'relative offset is {rel_offset}')
                # resize single frame and mask
                self.resize_frame_and_mask(iinfo, finalshape, window, scale)
                self.logger.debug('---')

    def resize_frame_and_mask(self, iinfo, finalshape, window, scale):
        self.logger.debug(f'resizing frame {iinfo.resized_base}')
        with iinfo.origin.open() as hdul:
            baseshape = hdul[0].data.shape

            # update CRPIX1 and CRPIX2 in header (to fix the WCS in the resized images)
            crpix1 = hdul[0].header['crpix1']
            crpix2 = hdul[0].header['crpix2']
            hdul[0].header['crpix1'] = crpix1 + iinfo.rel_offset[1]
            hdul[0].header['crpix2'] = crpix2 + iinfo.rel_offset[0]

            # resize_fits saves the resized image in frame_name
            frame_name = iinfo.resized_base
            resize_fits(hdul, frame_name, finalshape, iinfo.valid_region,
                        window=window, scale=scale, dtype='float32')

        self.logger.debug('resizing mask  %s', iinfo.resized_mask)
        if iinfo.mask is None:
            self.logger.warning('BPM missing, use zeros instead')
            false_mask = numpy.zeros(baseshape, dtype='int16')
            hdul_dum = fits.HDUList(fits.PrimaryHDU(false_mask))
            iinfo.mask = hdul_dum
        elif isinstance(iinfo.mask, nfcom.Extension):
            ename = iinfo.mask.name
            with iinfo.origin.open() as hdul:
                iinfo.mask = fits.HDUList(hdul[ename].copy())

        # We don't conserve the sum of the values of the frame here, just
        # expand the mask
        mask_name = iinfo.resized_mask
        resize_fits(iinfo.mask, mask_name, finalshape, iinfo.valid_region,
                    fill=1, window=window, scale=scale, conserve=False)

    def create_objmask(self, img, seeing_fwhm, step=0):

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

            # Take the number of combined images from the combined image (extension MAP)
            wm = img['MAP'].data.copy()
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

        self.logger.info('Running source extraction in previous result')
        objects, objmask = sep.extract(data_sub, 1.5, err=bkg.globalrms,
                                       mask=border, segmentation_map=True)
        self.logger.debug(f'... saving segmentation mask: {name_segmask(step)}')
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
                             maxsep_time=5.0,
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
                                  distance_upper_bound=maxsep_time * SCALE)

        nsky = len(sarray)

        for tid, idss in enumerate(idxs):
            self.logger.info('---')
            self.logger.info(f'image {tid + 1} / {len(idxs)}')
            tf = targetframes[tid]
            self.logger.info(f"Step {step}, SC: computing advanced sky for {tf.label} using '{method.__name__}'")
            try:
                # filter(lambda x: x < nsky, idss)
                locskyframes = []
                for si in idss:
                    if tid == si:
                        # this sky frame it is the current frame, reject
                        continue
                    if si < nsky:
                        self.logger.debug(f'Step {step}, SC: {skyframes[si].label} is a sky frame')
                        locskyframes.append(skyframes[si])
                self.compute_advanced_sky_for_frame(
                    tf, locskyframes, step=step,
                    method=method, method_kwargs=method_kwargs,
                    nside_adhoc_sky_correction=nside_adhoc_sky_correction
                )
            except IndexError:
                self.logger.error(f'No sky image available for frame {tf.lastname}')
                raise

    def compute_advanced_sky_for_frame(self, frame, skyframes, step=0,
                                       method=None, method_kwargs=None,
                                       nside_adhoc_sky_correction=0):
        self.logger.info('Correcting sky in frame %s', frame.lastname)
        self.logger.info('with sky computed from frames:')
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

                if i.objmask_data is not None:
                    msk = i.objmask_data
                    masks.append(msk)
                    self.logger.debug('object mask (including footprint) is shared')
                elif i.objmask is not None:
                    hdulistmask = fits.open(
                        i.objmask, mode='readonly', memmap=True)
                    # note that this image has the correct shape (2048x2048)
                    # and there is no need to use i.valid_region; in addition,
                    # this image also contain the footprint
                    msk = hdulistmask['primary'].data
                    masks.append(msk)
                    desc.append(hdulistmask)
                    self.logger.debug('object mask is particular (it must contain the footprint)')
                    self.logger.debug(f'reading {i.objmask}')
                else:
                    footprint = i.mask[0].data
                    msk = footprint
                    self.logger.warning(f'no object mask (using only footprint) for {filename}')

                if msk is None:   # this should never happen now (after including the footprint)
                    scales.append(numpy.median(data[-1]))
                else:
                    scales.append(numpy.median(data[-1][msk == 0]))

            self.logger.debug(f"computing scaled background with {len(data)} frames using '{method.__name__}'")
            self.logger.debug(f"... scales: {scales}")
            # note: this sky is scaled to have a mean value of 1.0 (using the unmasked pixels)
            sky, _, num = method(data, masks, scales=scales, **method_kwargs)

            with fits.open(frame.lastname) as hdulist:
                data = hdulist['primary'].data
                valid = data[frame.valid_region]

                if frame.objmask_data is not None:
                    self.logger.debug('object mask defined (including footprint)')
                    msk = frame.objmask_data
                else:
                    self.logger.debug('object mask empty (using only footprint)')
                    footprint = frame.mask[0].data
                    msk = footprint

                skymedian = numpy.median(valid[msk == 0])
                self.logger.debug(f'rescaling background with skymedian {skymedian}')

            # avoid pixels without sky information
            if numpy.any(num == 0):
                self.logger.warning('pixels without sky information found (set to skymedian)')
                sky[num == 0] = 1.0

            # rescale sky to have a mean value equal to skymedian
            sky *= skymedian

        finally:
            # Closing all FITS files
            for hdl in desc:
                hdl.close()

        # the following code is not necessary because we have already avoided
        # those pixels where num == 0
        """
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
            self.logger.debug(f'saving sky background mask {name}')
            fits.writeto(name, binmask.astype('int16'), overwrite=True)
        """

        name_sky = name_skybackground(frame.label, step)
        self.logger.debug(f'saving sky background {name_sky}')
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
            self.logger.debug(f'saving sky subtracted result {frame.lastname}')

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
