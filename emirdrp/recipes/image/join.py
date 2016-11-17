#
# Copyright 2014-2016 Universidad Complutense de Madrid
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

"""Dither Image recipe for EMIR"""

from __future__ import division, print_function


import datetime
import uuid

import numpy
from astropy.io import fits
from numina.core import Product, RecipeError
from numina.core.requirements import ObservationResultRequirement
from numina.array import combine
from numina.array import combine_shape, combine_shapes
from numina.array import resize_arrays, resize_arrays_alt
from numina.array.utils import coor_to_pix, image_box2d
from numina.core import ObservationResult
from numina.flow.processing import SkyCorrector

from emirdrp.processing.wcs import offsets_from_wcs, reference_pix_from_wcs
from emirdrp.processing.corr import offsets_from_crosscor
from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.ext.gtc import RUN_IN_GTC
from emirdrp.processing.combine import segmentation_combined
import emirdrp.decorators


class JoinDitheredImagesRecipe(EmirRecipe):
    """Combine single exposures obtained in dithered mode"""

    obresult = ObservationResultRequirement()
    frame = Product(DataFrameType)
    sky = Product(DataFrameType, optional=True)
    #
    # Accumulate Frame results
    accum = Product(DataFrameType, optional=True)

    @classmethod
    def build_recipe_input(cls, obsres, dal, pipeline='default'):
        if RUN_IN_GTC:
            cls.logger.debug('Using GTC version of build_recipe_input in DitheredImages')
            return cls.build_recipe_input_gtc(obsres, dal, pipeline=pipeline)
        else:
            return super(JoinDitheredImagesRecipe, cls).build_recipe_input(obsres, dal, pipeline=pipeline)

    @classmethod
    def build_recipe_input_gtc(cls, obsres, dal, pipeline='default'):

        # FIXME: this method will work only in GTC
        # stareImagesIds = obsres['stareImagesIds']._v
        stareImagesIds = obsres.stareImagesIds
        cls.logger.info('obsres(%s): %s', type(obsres), dir(obsres))
        cls.logger.info('STARE IMAGES IDS: %s', stareImagesIds)
        stareImages = []
        # Field to query the results
        key_field = 'frame'
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            # This 'frame' is the name of the product in RecipeResult
            # there is also a 'sky' field
            elements = subres['elements']
            stareImages.append(elements[key_field])

        newOR = ObservationResult()
        newOR.frames = stareImages
        # obsres['obresult'] = newOR
        # print('Adding RI parameters ', obsres)
        # newRI = DitheredImageARecipeInput(**obsres)
        newRI = cls.create_input(obresult=newOR)
        return newRI

    #@emirdrp.decorators.aggregate
    @emirdrp.decorators.loginfo
    def run(self, rinput):
        partial_result = self.run_single(rinput)
        new_result = self.aggregate_result(partial_result, rinput)
        return new_result

    def aggregate_result(self, partial_result, rinput):

        obresult = rinput.obresult
        # Check if this is our first run
        naccum = getattr(obresult, 'naccum', 0)
        accum = getattr(obresult, 'accum', None)

        frame = partial_result.frame

        if naccum == 0:
            self.logger.debug('naccum is not set, do not accumulate')
            return partial_result
        elif naccum == 1:
            self.logger.debug('round %d initialize accumulator', naccum)
            newaccum = frame
        elif naccum > 1:
            self.logger.debug('round %d of accumulation', naccum)
            newaccum = self.aggregate_frames(accum, frame, naccum)
        else:
            msg = 'naccum set to %d, invalid' % (naccum, )
            self.logger.error(msg)
            raise RecipeError(msg)

        # Update partial result
        partial_result.accum = newaccum

        return partial_result

    def aggregate_frames(self, accum, frame, naccum):
        return self.aggregate2(accum, frame, naccum)

    def run_single(self, rinput):

        # Open all images
        obresult = rinput.obresult

        data_hdul = []
        for f in obresult.frames:
            img = f.open()
            data_hdul.append(img)

        use_errors = True
        # Initial checks
        baseimg = data_hdul[0]
        has_num_ext = 'NUM' in baseimg
        has_bpm_ext = 'BPM' in baseimg
        baseshape = baseimg[0].shape
        subpixshape = baseshape
        base_header = baseimg[0].header
        compute_sky = 'NUM-SK' not in base_header

        self.logger.debug('base image is: %s', self.datamodel.get_imgid(img))
        self.logger.debug('images have NUM extension: %s', has_num_ext)
        self.logger.debug('images have BPM extension: %s', has_bpm_ext)
        self.logger.debug('compute sky is needed: %s', compute_sky)

        self.logger.info('Computing offsets from WCS information')

        finalshape, offsetsp, refpix = self.compute_offset_wcs(
            rinput.obresult.frames,
            baseshape,
            subpixshape
        )

        self.logger.debug("Relative offsetsp %s", offsetsp)
        self.logger.info('Shape of resized array is %s', finalshape)

        # Resizing target frames
        data_arr_r, regions = resize_arrays(
            [m[0].data for m in data_hdul],
            subpixshape,
            offsetsp,
            finalshape,
            fill=1
        )

        try:
            self.logger.debug("Compute cross-correlation of images")
            # A square of 100x100 in the center of the image
            xref_cross = finalshape[1] // 2
            yref_cross = finalshape[0] // 2
            box = 50
            self.logger.debug("Reference position is %d  %d", xref_cross + 1, yref_cross + 1)
            self.logger.debug("Reference regions is %d", 2 * box + 1)
            region = image_box2d(xref_cross, yref_cross, finalshape, (box, box))
            finalshape2, offsetsp2 = self.compute_offset_crosscor(data_arr_r, region, finalshape)
            self.logger.debug("Relative offsetsp (crosscorr) %s", offsetsp2)
            self.logger.info('Shape of resized array (crosscorr) is %s', finalshape2)
        except Exception as error:
            self.logger.warning('Error during cross-correlation, %s', error)

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
            masks = [false_mask for _ in data_arr_r]

        self.logger.debug('resize bad pixel masks')
        mask_arr_r, _ = resize_arrays(masks, subpixshape, offsetsp, finalshape, fill=1)

        if compute_sky:
            self.logger.debug("compute sky")

            omasks = self.compute_object_masks(
                data_arr_r,
                mask_arr_r,
                has_bpm_ext,
                regions,
                masks
            )

            sky_result = self.compute_sky(data_hdul, omasks, base_header, use_errors)
            sky_data = sky_result[0].data
            self.logger.debug('sky image has shape %s', sky_data.shape)

            self.logger.info('sky correction in individual images')
            corrector = SkyCorrector(
                sky_data,
                self.datamodel,
                calibid=self.datamodel.get_imgid(sky_result)
            )
            # If we do not update keyword SKYADD
            # there is no sky subtraction
            for m in data_hdul:
                m[0].header['SKYADD'] = True
            # this is a little hackish
            data_hdul_s = [corrector(m) for m in data_hdul]
            #data_arr_s = [m[0].data - sky_data for m in data_hdul]
            base_header = data_hdul_s[0][0].header
            self.logger.info('resize sky-corrected images')
            data_arr_sr, _ = resize_arrays(
                [f[0].data for f in data_hdul_s],
                subpixshape,
                offsetsp,
                finalshape,
                fill=0
            )
        else:
            self.logger.debug("not computing sky")
            sky_result = None
            data_arr_sr = data_arr_r

        # Position of refpixel in final image
        refpix_final = refpix + offsetsp[0]
        self.logger.info('Position of refpixel in final image %s', refpix_final)

        self.logger.info('Combine target images (final)')
        method = combine.median
        out = method(data_arr_sr, masks=mask_arr_r, dtype='float32')

        self.logger.debug('create result image')
        hdu = fits.PrimaryHDU(out[0], header=base_header)
        self.logger.debug('update result header')
        hdr = hdu.header
        self.set_base_headers(hdr)
        hdr['IMGOBBL'] = 0
        hdr['TSUTC2'] = data_hdul[-1][0].header['TSUTC2']
        # Update obsmode in header
        hdr['OBSMODE'] = 'DITHERED_IMAGE'
        hdu.header['history'] = "Combined %d images using '%s'" % (
            len(data_hdul),
            method.__name__
        )
        hdu.header['history'] = 'Combination time {}'.format(
            datetime.datetime.utcnow().isoformat()
        )
        # Update NUM-NCOM, sum of individual frames
        ncom = 0
        for hdul in data_hdul:
            ncom += hdul[0].header.get('NUM-NCOM', 1)
        hdr['NUM-NCOM'] = ncom
        # Update WCS, approximate solution
        hdr['CRPIX1'] += offsetsp[0][0]
        hdr['CRPIX2'] += offsetsp[0][1]

        #
        if use_errors:
            varhdu = fits.ImageHDU(out[1], name='VARIANCE')
            num = fits.ImageHDU(out[2], name='MAP')
            hdulist = fits.HDUList([hdu, varhdu, num])
        else:
            hdulist = fits.HDUList([hdu])

        result = self.create_result(frame=hdulist, sky=sky_result)
        self.logger.info('end of dither recipe')
        return result

    def compute_offset_wcs(self, frames, baseshape, subpixshape):

        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2).astype('float')
        offsets_xy = offsets_from_wcs(frames, refpix)
        self.logger.debug("offsets_xy %s", offsets_xy)
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:, ::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        self.logger.info('Computing relative offsets')
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)

        return finalshape, offsetsp, refpix

    def compute_offset_crosscor(self, arrs, region, subpixshape):
        offsets_xy = offsets_from_crosscor(arrs, region, refine=True, order='xy')
        self.logger.debug("offsets_xy cross-corr %s", offsets_xy)
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:, ::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        self.logger.info('Computing relative offsets from cross-corr')
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)

        return finalshape, offsetsp

    def compute_shapes_wcs(self, frames):

        # Better near the center...
        shapes = [frame.open()[0].shape for frame in frames]
        ref_pix_xy_0 = (shapes[0][1] // 2, shapes[0][0] // 2)
        #
        ref_coor_xy = reference_pix_from_wcs(frames, ref_pix_xy_0)
        # offsets_xy = offsets_from_wcs(frames, numpy.asarray([ref_pix_xy_0]))
        # ll = [(-a[0]+ref_coor_xy[0][0], -a[1]+ref_coor_xy[0][1]) for a in ref_coor_xy]

        self.logger.debug("ref_coor_xy %s", ref_coor_xy)
        # Transform to pixels, integers
        ref_pix_xy = [coor_to_pix(c, order='xy') for c in ref_coor_xy]

        self.logger.info('Computing relative shapes')
        finalshape, partialshapes, finalpix_xy = combine_shapes(shapes, ref_pix_xy, order='xy')

        return finalshape, partialshapes, ref_pix_xy_0, finalpix_xy

    def compute_object_masks(self, data_arr_r, mask_arr_r, has_bpm_ext, regions, masks):

        method = combine.mean

        self.logger.info(
            "initial stacking, %d images, with offsets using '%s'",
            len(data_arr_r),
            method.__name__
        )
        data1 = method(data_arr_r, masks=mask_arr_r, dtype='float32')

        self.logger.info('obtain segmentation mask')
        segmap = segmentation_combined(data1[0])
        # submasks
        if not has_bpm_ext:
            omasks = [(segmap[region] > 0) for region in regions]
        else:
            omasks = [((segmap[region] > 0) & bpm) for region, bpm in zip(regions, masks)]

        return omasks

    def compute_sky(self, data_hdul, omasks, base_header, use_errors):
        method = combine.mean

        self.logger.info('recombine images with segmentation mask')
        sky_data = method([m[0].data for m in data_hdul], masks=omasks, dtype='float32')

        hdu = fits.PrimaryHDU(sky_data[0], header=base_header)
        points_no_data = (sky_data[2] == 0).sum()

        self.logger.debug('update created sky image result header')
        skyid = uuid.uuid1().hex
        hdu.header['UUID'] = skyid
        hdu.header['history'] = "Combined {} images using '{}'".format(
            len(data_hdul),
            method.__name__
        )

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

    def aggregate2(self, img1, img2, naccum):
        # FIXME, this is almost identical to run_single
        frames = [img1, img2]
        use_errors = True
        # Initial checks
        fframe = frames[0]
        img = fframe.open()
        base_header = img[0].header

        data_hdul = []
        for f in frames:
            img = f.open()
            data_hdul.append(img)

        self.logger.info('Computing offsets from WCS information')

        finalshape, partial_shapes, refpix_xy_0, refpix_final_xy = self.compute_shapes_wcs(frames)

        self.logger.info('Shape of resized array is %s', finalshape)
        self.logger.debug("partial shapes %s", partial_shapes)

        masks = []
        self.logger.debug('Obtains masks')
        for m in data_hdul:
            if 'NUM' in m:
                self.logger.debug('Using NUM extension as mask')
                mask = numpy.where(m['NUM'].data, 0, 1).astype('int16')
            elif 'BPM' in m:
                self.logger.debug('Using BPM extension as mask')
                mask = numpy.where(m['BPM'].data, 1, 0).astype('int16')
            else:
                self.logger.warning('BPM missing, use zeros instead')
                mask = numpy.zeros_like(m[0].data)
            masks.append(mask)

        # Resizing target frames
        data_arr_r = resize_arrays_alt(
            [m[0].data for m in data_hdul],
            partial_shapes,
            finalshape,
            fill=1
        )

        self.logger.debug('resize bad pixel masks')
        mask_arr_r = resize_arrays_alt(masks, partial_shapes, finalshape, fill=1)

        self.logger.debug("not computing sky")
        data_arr_sr = data_arr_r

        self.logger.info('Combine target images (final, aggregate)')
        self.logger.debug("weights for 'accum' and 'frame'")

        weight_accum = 2 * (1 - 1.0 / naccum)
        weight_frame = 2.0 / naccum
        scales = [1.0 / weight_accum, 1.0 / weight_frame]
        self.logger.debug("weights for 'accum' and 'frame', %s", scales)
        method = combine.mean

        out = method(data_arr_sr, masks=mask_arr_r, scales=scales, dtype='float32')

        self.logger.debug('create result image')
        hdu = fits.PrimaryHDU(out[0], header=base_header)
        self.logger.debug('update result header')
        hdr = hdu.header
        self.set_base_headers(hdr)
        hdr['IMGOBBL'] = 0
        hdr['TSUTC2'] = data_hdul[-1][0].header['TSUTC2']
        # Update obsmode in header
        hdr['OBSMODE'] = 'DITHERED_IMAGE'
        hdu.header['history'] = "Combined %d images using '%s'" % (
            len(data_hdul),
            method.__name__
        )
        hdu.header['history'] = 'Combination time {}'.format(
            datetime.datetime.utcnow().isoformat()
        )
        # Update NUM-NCOM, sum of individual frames
        ncom = 0
        for hdul in data_hdul:
            ncom += hdul[0].header['NUM-NCOM']
        hdr['NUM-NCOM'] = ncom
        # Update WCS, approximate solution
        hdr['CRPIX1'] += (refpix_final_xy[0] - refpix_xy_0[0])
        hdr['CRPIX2'] += (refpix_final_xy[1] - refpix_xy_0[1])

        #
        if use_errors:
            varhdu = fits.ImageHDU(out[1], name='VARIANCE')
            num = fits.ImageHDU(out[2], name='MAP')
            hdulist = fits.HDUList([hdu, varhdu, num])
        else:
            hdulist = fits.HDUList([hdu])

        return hdulist

    def offsets_from_crosscorr(self, arrs, region):
        ref_arr = arrs[0]
        arr1 = ref_arr[region]
        for idx, arr in enumerate(arrs[1:]):
            self.logger.debug('Cross-correlate image idx %d', idx)
            arr2 = arr[region]
