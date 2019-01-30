#
# Copyright 2008-2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Recipe for the reduction of imaging mode observations."""

import os
import shutil
import datetime
import uuid

import numpy
import sep
from astropy.io import fits
from scipy.spatial import KDTree as KDTree

from numina.core import Result, Requirement, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.core.query import ResultOf
import numina.array as narray
import numina.array.utils as nautils
import numina.array.combine as nacom
import numina.frame.combine as nfcom
import numina.processing as proc
from numina.frame import resize_fits, custom_region_to_str

import emirdrp.requirements as reqs
import emirdrp.products as prods
from emirdrp.instrument.channels import FULL
from emirdrp.processing.wcs import offsets_from_wcs_imgs
from emirdrp.processing.corr import offsets_from_crosscor, offsets_from_crosscor_regions
from emirdrp.core.recipe import EmirRecipe

from .naming import (name_redimensioned_frames, name_object_mask,
                     name_skybackground)
from .naming import name_skybackgroundmask, name_skysub_proc, name_skyflat
from .naming import name_skyflat_proc, name_segmask


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
    obresult = ObservationResultRequirement(query_opts=ResultOf('result_image', node='children'))

    master_bpm = reqs.MasterBadPixelMaskRequirement()

    offsets = Requirement(
        prods.CoordinateList2DType,
        'List of pairs of offsets',
        optional=True
    )
    refine_offsets = Parameter(False, 'Refine offsets by cross-correlation')
    iterations = Parameter(2, 'Iterations of the recipe')
    extinction = Parameter(0.0, 'Mean atmospheric extinction')

    sky_images = Parameter(
        5, 'Images used to estimate the '
           'background before and after current image')

    sky_images_sep_time = Parameter(
        10, 'Maximum time interval between target and sky images [minutes]'
    )

    result_image = Result(prods.ProcessedImage)
    result_sky = Result(prods.ProcessedImage, optional=True)

    def run(self, rinput):

        target_is_sky = True
        obresult = rinput.obresult
        sky_images = rinput.sky_images
        sky_images_sep_time = rinput.sky_images_sep_time
        baseshape = (2048, 2048)
        user_offsets = rinput.offsets
        extinction = rinput.extinction

        images_info = self.initial_classification(obresult, target_is_sky)

        # Resizing target frames
        target_info = [iinfo for iinfo in images_info if iinfo.valid_target]
        finalshape, offsetsp, refpix, offset_fc0 = self.compute_size(
            target_info, baseshape, user_offsets
        )

        self.resize(target_info, baseshape, offsetsp, finalshape)

        result = self.process_basic(images_info, target_is_sky=target_is_sky,
                                     extinction=extinction)

        if rinput.refine_offsets:
            self.logger.debug("Compute cross-correlation of images")
            # regions_c = self.compute_regions(finalshape, box=200, corners=True)

            # Regions frm bright objects
            regions_c = self.compute_regions_from_objs(result[0].data, finalshape, box=40)

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
                self.logger.debug('Total offsets: %s', offsets_xy_t)
                self.logger.info('Computing relative offsets from cross-corr')
                finalshape2, offsetsp2 = narray.combine_shape(baseshape, offsets_fc_t)
                #
                self.logger.debug("Relative offsetsp (crosscorr) %s", offsetsp2)
                self.logger.info('Shape of resized array (crosscorr) is %s', finalshape2)

                # Resizing target imgs
                self.logger.debug("Resize to final offsets")
                self.resize(target_info, baseshape, offsetsp2, finalshape2)
            except Exception as error:
                self.logger.warning('Error during cross-correlation, %s', error)

        result = self.process_basic(images_info, target_is_sky=target_is_sky,
                                     extinction=extinction)

        step = 1

        while step <= rinput.iterations:
            result = self.process_advanced(
                images_info, result, step, target_is_sky,
                maxsep=sky_images_sep_time, nframes=sky_images,
                extinction=extinction
            )
            step += 1

        return self.create_result(result_image=result)

    def compute_offset_xy_crosscor_regions(self, iinfo, regions, refine=False, tol=0.5):

        names = [frame.lastname for frame in iinfo]
        print(names)
        print(regions)
        with nfcom.manage_fits(names) as imgs:
            arrs = [img[0].data for img in imgs]
            offsets_xy = offsets_from_crosscor_regions(
                arrs, regions,
                refine=refine, order='xy', tol=tol
            )
            self.logger.debug("offsets_xy cross-corr %s", offsets_xy)
            # Offsets in numpy order, swaping
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
            self.logger.debug('Computing offsets from WCS information')
            with nfcom.manage_fits(img.origin for img in target_info) as images:
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
        self.logger.debug("Relative offsetsp %s", offsetsp)
        self.logger.info('Shape of resized array is %s', finalshape)
        return finalshape, offsetsp, refpix, list_of_offsets


    def process_basic(self, images_info, target_is_sky=True, extinction=0.0):

        step = 0

        target_info = [iinfo for iinfo in images_info if iinfo.valid_target]
        sky_info = [iinfo for iinfo in images_info if iinfo.valid_sky]

        self.logger.info("Step %d, SF: compute superflat", step)
        sf_arr = self.compute_superflat(images_info)

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
        result = self.combine_frames(target_info, extinction=extinction)
        self.logger.info('Step %d, finished', step)

        return result


    def process_advanced(self, images_info, result, step, target_is_sky=True,
                maxsep=5.0, nframes=6, extinction=0):

        seeing_fwhm = None
        baseshape = (2048, 2048)
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
        sf_arr = self.compute_superflat(
            sky_info, segmask=objmask, step=step)

        # Apply superflat
        self.logger.info("Step %d, SF: apply superflat", step)
        for iinfo in images_info:
            self.correct_superflat(iinfo, sf_arr, step=step, save=True)

        self.logger.info('Step %d, advanced sky correction (SC)', step)
        self.compute_advanced_sky(target_info, objmask,
                                  skyframes=sky_info,
                                  target_is_sky=target_is_sky,
                                  maxsep=maxsep,
                                  nframes=nframes,
                                  step=step)

        # Combining the images
        self.logger.info("Step %d, Combining the images", step)
        # FIXME: only for science
        result = self.combine_frames(
            target_info, extinction, step=step)
        return result


    def compute_simple_sky_for_frame(self, frame, skyframe, step=0, save=True):
        self.logger.info('Correcting sky in frame %s', frame.lastname)
        self.logger.info('with sky computed from frame %s', skyframe.lastname)

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
            self.logger.debug('images have BPM extension: %s', has_bpm_ext)

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
                    #targetframes.append(iinfo)
                    if target_is_sky:
                        iinfo.valid_sky = True
                        #skyframes.append(iinfo)
                if iinfo.itype == 'SKY':
                    iinfo.valid_sky = True
                    #skyframes.append(iinfo)
                images_info.append(iinfo)

        return images_info


    def compute_superflat(self, images_info, segmask=None, step=0):

        self.logger.info("Step %d, SF: combining the frames without offsets", step)

        base_imgs = [img.resized_base for img in images_info]
        with nfcom.manage_fits(base_imgs) as imgs:

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
                    self.logger.debug('Step %d, opening resized mask %s',
                                      step, frame.resized_mask)
                    hdulist = fits.open(
                         frame.resized_mask, memmap=True, mode='readonly')
                    #filelist.append(hdulist)
                    masks.append(hdulist['primary'].data[frame.valid_region])
                masks = None

            self.logger.debug('Step %d, combining %d frames', step, len(data))
            sf_data, _sf_var, sf_num = nacom.median(data, masks, scales=scales,
                                                 dtype='float32',
                                                 #blank=1.0 / scales[0]
                                                 )

        # Normalize, flat has mean = 1
        sf_data[sf_data == 0] = 1e-5
        sf_data /= sf_data.mean()
        #sf_data[sf_data <= 0] = 1.0

        # Auxiliary data
        sfhdu = fits.PrimaryHDU(sf_data)
        self.save_intermediate_img(sfhdu, name_skyflat('comb', step))
        return sf_data


    def run_single(self, rinput):

        # FIXME: remove this, is deprecated

        obresult = rinput.obresult

        # just in case images are in result, instead of frames
        if not obresult.frames:
            frames = obresult.results
        else:
            frames = obresult.frames

        img_info = []
        data_hdul = []
        for f in frames:
            img = f.open()
            data_hdul.append(img)
            info = {}
            info['tstamp'] = img[0].header['tstamp']
            info['airmass'] = img[0].header['airmass']
            img_info.append(info)

        channels = FULL

        use_errors = True
        # Initial checks
        baseimg = data_hdul[0]
        has_num_ext = 'NUM' in baseimg
        has_bpm_ext = 'BPM' in baseimg
        baseshape = baseimg[0].shape
        subpixshape = baseshape
        base_header = baseimg[0].header
        compute_sky = 'NUM-SK' not in base_header
        compute_sky_advanced = False

        self.logger.debug('base image is: %s', self.datamodel.get_imgid(baseimg))
        self.logger.debug('images have NUM extension: %s', has_num_ext)
        self.logger.debug('images have BPM extension: %s', has_bpm_ext)
        self.logger.debug('compute sky is needed: %s', compute_sky)

        if compute_sky:
            self.logger.info('compute sky simple')
            sky_result = self.compute_sky_simple(data_hdul, use_errors=False)
            self.save_intermediate_img(sky_result, 'sky_init.fits')
            sky_result.writeto('sky_init.fits', overwrite=True)
            sky_data = sky_result[0].data
            self.logger.debug('sky image has shape %s', sky_data.shape)

            self.logger.info('sky correction in individual images')
            corrector = proc.SkyCorrector(
                sky_data,
                self.datamodel,
                calibid=self.datamodel.get_imgid(sky_result)
            )
            # If we do not update keyword SKYADD
            # there is no sky subtraction
            for m in data_hdul:
                m[0].header['SKYADD'] = True
            # this is a little hackish
            # sky corrected
            data_hdul_s = [corrector(m) for m in data_hdul]
            base_header = data_hdul_s[0][0].header
        else:
            sky_result = None
            data_hdul_s = data_hdul

        self.logger.info('Computing offsets from WCS information')

        finalshape, offsetsp, refpix, offset_xy0 = self.compute_offset_wcs_imgs(
            data_hdul_s,
            baseshape,
            subpixshape
        )

        self.logger.debug("Relative offsetsp %s", offsetsp)
        self.logger.info('Shape of resized array is %s', finalshape)

        # Resizing target imgs
        data_arr_sr, regions = narray.resize_arrays(
            [m[0].data for m in data_hdul_s],
            subpixshape,
            offsetsp,
            finalshape,
            fill=1
        )

        if has_num_ext:
            self.logger.debug('Using NUM extension')
            masks = [numpy.where(m['NUM'].data, 0, 1).astype('int16') for m in data_hdul]
        elif has_bpm_ext:
            self.logger.debug('Using BPM extension')
            #
            masks = [numpy.where(m['BPM'].data, 1, 0).astype('int16') for m in data_hdul]
        else:
            self.logger.warning('BPM missing, use zeros instead')
            false_mask = numpy.zeros(baseshape, dtype='int16')
            masks = [false_mask for _ in data_arr_sr]

        self.logger.debug('resize bad pixel masks')
        mask_arr_r, _ = narray.resize_arrays(masks, subpixshape, offsetsp, finalshape, fill=1)

        if self.intermediate_results:
            self.logger.debug('save resized intermediate img')
            for idx, arr_r in enumerate(data_arr_sr):
                self.save_intermediate_array(arr_r, 'interm1_%03d.fits' % idx)

        hdulist = self.combine2(data_arr_sr, mask_arr_r, data_hdul, offsetsp, use_errors)

        self.save_intermediate_img(hdulist, 'result_initial1.fits')

        compute_cross_offsets = True
        if compute_cross_offsets:

            self.logger.debug("Compute cross-correlation of images")
            # regions_c = self.compute_regions(finalshape, box=200, corners=True)

            # Regions frm bright objects
            regions_c = self.compute_regions_from_objs(hdulist[0].data, finalshape, box=20)

            try:

                offsets_xy_c = self.compute_offset_xy_crosscor_regions(
                    data_arr_sr, regions_c, refine=True, tol=1
                )
                #
                # Combined offsets
                # Offsets in numpy order, swaping
                offsets_xy_t = offset_xy0 - offsets_xy_c
                offsets_fc = offsets_xy_t[:, ::-1]
                offsets_fc_t = numpy.round(offsets_fc).astype('int')
                self.logger.debug('Total offsets: %s', offsets_xy_t)
                self.logger.info('Computing relative offsets from cross-corr')
                finalshape, offsetsp = narray.combine_shape(subpixshape, offsets_fc_t)
                #
                self.logger.debug("Relative offsetsp (crosscorr) %s", offsetsp)
                self.logger.info('Shape of resized array (crosscorr) is %s', finalshape)

                # Resizing target imgs
                self.logger.debug("Resize to final offsets")
                data_arr_sr, regions = narray.resize_arrays(
                    [m[0].data for m in data_hdul_s],
                    subpixshape,
                    offsetsp,
                    finalshape,
                    fill=1
                )

                if self.intermediate_results:
                    self.logger.debug('save resized intermediate2 img')
                    for idx, arr_r in enumerate(data_arr_sr):
                        self.save_intermediate_array(arr_r, 'interm2_%03d.fits' % idx)

                self.logger.debug('resize bad pixel masks')
                mask_arr_r, _ = narray.resize_arrays(masks, subpixshape, offsetsp, finalshape, fill=1)

                hdulist = self.combine2(data_arr_sr, mask_arr_r, data_hdul, offsetsp, use_errors)

                self.save_intermediate_img(hdulist, 'result_initial2.fits')
            except Exception as error:
                self.logger.warning('Error during cross-correlation, %s', error)


        catalog, objmask = self.create_object_catalog(hdulist[0].data, border=50)

        data_arr_sky = [sky_result[0].data for _ in data_arr_sr]
        data_arr_0 = [(d[r] + s) for d, r, s in zip(data_arr_sr, regions, data_arr_sky)]
        data_arr_r = [d.copy() for d in data_arr_sr]

        for inum in range(1, rinput.iterations + 1):
            # superflat
            sf_data = self.compute_superflat(data_arr_0, objmask, regions, channels)
            fits.writeto('superflat_%d.fits' % inum, sf_data, overwrite=True)
            # apply superflat
            data_arr_rf = data_arr_r
            for base, arr, reg in zip(data_arr_rf, data_arr_0, regions):
                arr_f = arr / sf_data
                #arr_f = arr
                base[reg] = arr_f

            # compute sky advanced
            data_arr_sky = []
            data_arr_rfs = []
            self.logger.info('Step %d, SC: computing advanced sky', inum)
            scale = rinput.sky_images_sep_time * 60
            tstamps = numpy.array([info['tstamp'] for info in img_info])
            for idx, hdu in enumerate(data_hdul):
                diff1 = tstamps - tstamps[idx]
                idxs1 = (diff1 > 0) & (diff1 < scale)
                idxs2 = (diff1 < 0) & (diff1 > -scale)
                l1, = numpy.nonzero(idxs1)
                l2, = numpy.nonzero(idxs2)
                limit1 = l1[-rinput.sky_images:]
                limit2 = l2[:rinput.sky_images]
                len_l1 =len(limit1)
                len_l2 = len(limit2)
                self.logger.info('For image %s, using %d-%d images)', idx,
                                 len_l1, len_l2)
                if len_l1 + len_l2 == 0:
                    self.logger.error(
                        'No sky image available for frame %d', idx)
                    raise ValueError('No sky image')
                skydata = []
                skymasks = []
                skyscales = []
                my_region = regions[idx]
                my_sky_scale = numpy.median(data_arr_rf[idx][my_region])
                for i in numpy.concatenate((limit1, limit2)):
                    region_s = regions[i]
                    data_s = data_arr_rf[i][region_s]
                    mask_s = objmask[region_s]
                    scale_s = numpy.median(data_s)
                    skydata.append(data_s)
                    skymasks.append(mask_s)
                    skyscales.append(scale_s)
                self.logger.debug('computing background with %d frames', len(skydata))
                sky, _, num = nacom.median(skydata, skymasks, scales=skyscales)
                # rescale
                sky *= my_sky_scale

                binmask = num == 0

                if numpy.any(binmask):
                    # We have pixels without
                    # sky background information
                    self.logger.warn('pixels without sky information when correcting %d',
                                 idx)

                    # FIXME: during development, this is faster
                    # sky[binmask] = sky[num != 0].mean()
                    # To continue we interpolate over the patches
                    narray.fixpix2(sky, binmask, out=sky, iterations=1)

                name = 'sky_%d_%03d.fits' % (inum, idx)
                fits.writeto(name, sky, overwrite=True)
                name = 'sky_binmask_%d_%03d.fits' % (inum, idx)
                fits.writeto(name, binmask.astype('int16'), overwrite=True)

                data_arr_sky.append(sky)
                arr = numpy.copy(data_arr_rf[idx])
                arr[my_region] = data_arr_rf[idx][my_region] - sky
                data_arr_rfs.append(arr)
                # subtract sky advanced

            if self.intermediate_results:
                self.logger.debug('save resized intermediate img')
                for idx, arr_r in enumerate(data_arr_rfs):
                    self.save_intermediate_array(arr_r, 'interm_%d_%03d.fits' % (inum, idx))


            hdulist = self.combine2(data_arr_rfs, mask_arr_r, data_hdul, offsetsp, use_errors)

            self.save_intermediate_img(hdulist, 'result_%d.fits' % inum)

            # For next step
            catalog, objmask = self.create_object_catalog(hdulist[0].data, border=50)

            data_arr_0 = [(d[r] + s) for d, r, s in zip(data_arr_rfs, regions, data_arr_sky)]
            data_arr_r = [d.copy() for d in data_arr_rfs]

        result = self.create_result(frame=hdulist)
        self.logger.info('end of dither recipe')
        return result

    def compute_sky_advanced(self, data_hdul, omasks, base_header, use_errors):
        method = narray.combine.mean

        self.logger.info('recombine images with segmentation mask')
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

    def combine_frames(self, frames, extinction, out=None, step=0):
        self.logger.debug('Step %d, opening sky-subtracted frames', step)

        def fits_open(name):
            """Open FITS with memmap in readonly mode"""
            return fits.open(name, mode='readonly', memmap=True)

        frameslll = [fits_open(frame.lastname)
                     for frame in frames if frame.valid_target]
        self.logger.debug('Step %d, opening mask frames', step)
        mskslll = [fits_open(frame.resized_mask)
                   for frame in frames if frame.valid_target]

        self.logger.debug('Step %d, combining %d frames', step, len(frameslll))
        try:
            extinc = [pow(10, -0.4 * frame.metadata['airmass'] * extinction)
                      for frame in frames if frame.valid_target]
            data = [i['primary'].data for i in frameslll]
            masks = [i['primary'].data for i in mskslll]
            headers = [i['primary'].header for i in frameslll]

            out = nacom.median(data, masks, scales=extinc, dtype='float32', out=out)

            base_header = headers[0]
            hdu = fits.PrimaryHDU(out[0], header=base_header)
            hdu.header['history'] = "Combined %d images using '%s'" % (len(frameslll), 'median')
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
                self.logger.debug('%s, valid region is %s, relative offset is %s',
                              frame.label, custom_region_to_str(region),
                              rel_offset
                              )
                self.resize_frame_and_mask(
                    frame, finalshape, framen, maskn, window, scale)

    def resize_frame_and_mask(self, frame, finalshape,
                              framen, maskn, window, scale):
        self.logger.info('Resizing frame %s', frame.label)
        with frame.origin.open() as hdul:
            baseshape = hdul[0].data.shape

            # FIXME: Resize_fits saves the resized image in framen
            resize_fits(hdul, framen, finalshape, frame.valid_region,
                        window=window, scale=scale, dtype='float32')

        self.logger.info('Resizing mask %s', frame.label)
        # We don't conserve the sum of the values of the frame here, just
        # expand the mask

        if frame.mask is None:
            self.logger.warning('BPM missing, use zeros instead')
            false_mask = numpy.zeros(baseshape, dtype='int16')
            hdum = fits.HDUList(fits.PrimaryHDU(false_mask))
            frame.mask = hdum #DataFrame(frame=hdum)
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
                             step=0, save=True):

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
                self.logger.info('Step %d, SC: computing advanced sky for %s',
                             step, tf.label)
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
                    tf, locskyframes, step=step, save=save)
            except IndexError:
                self.logger.error(
                    'No sky image available for frame %s', tf.lastname)
                raise

    def compute_advanced_sky_for_frame(self, frame, skyframes,
                                       step=0, save=True):
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
                #scales.append(numpy.median(data[-1]))
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

            self.logger.debug('computing background with %d frames', len(data))
            sky, _, num = nacom.median(data, masks)#, scales=scales)

        finally:
            # Closing all FITS files
            for hdl in desc:
                hdl.close()

        if numpy.any(num == 0):
            # We have pixels without
            # sky background information
            self.logger.warn('pixels without sky information when correcting %s',
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


    def compute_regions_from_objs(self, arr, finalshape, box=50, corners=True):
        regions = []
        catalog, mask = self.create_object_catalog(arr, border=300)

        self.save_intermediate_array(mask, 'objmask.fits')
        # with the catalog, compute 5 objects

        LIMIT_AREA = 5000
        NKEEP = 1
        idx_small = catalog['npix'] < LIMIT_AREA
        objects_small = catalog[idx_small]
        idx_flux = objects_small['flux'].argsort()
        objects_nth = objects_small[idx_flux][-NKEEP:]
        for obj in objects_nth:
            print('ref is', obj['x'], obj['y'])
            region = nautils.image_box2d(obj['x'], obj['y'], finalshape, (box, box))
            print(region)
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
