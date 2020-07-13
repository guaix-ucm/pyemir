#
# Copyright 2014-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""Dither Image recipe for EMIR"""

from __future__ import division, print_function


import datetime
import uuid

import numpy
import sep
from astropy.io import fits
from numina.core import Result, Requirement, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.array import combine
from numina.array import combine_shape, combine_shapes
from numina.array import resize_arrays, resize_arrays_alt
from numina.array.combine import flatcombine, median, quantileclip
from numina.array.utils import coor_to_pix, image_box2d
import numina.processing as proc
from numina.core.query import ResultOf
from numina.array import fixpix2
from numina.frame.utils import copy_img

from emirdrp.instrument.channels import FULL
import emirdrp.products as prods
import emirdrp.requirements as reqs
import emirdrp.decorators
from emirdrp.processing.wcs import offsets_from_wcs_imgs, reference_pix_from_wcs_imgs
from emirdrp.processing.corr import offsets_from_crosscor, offsets_from_crosscor_regions
from emirdrp.core.recipe import EmirRecipe
from emirdrp.processing.combine import segmentation_combined


class JoinDitheredImagesRecipe(EmirRecipe):
    """Combine single exposures obtained in dithered mode"""

    obresult = ObservationResultRequirement(query_opts=ResultOf(
        'STARE_IMAGE.reduced_image', node='children', id_field="stareImagesIds"))
    accum_in = Requirement(prods.ProcessedImage,
                           description='Accumulated result',
                           optional=True,
                           destination='accum',
                           query_opts=ResultOf('DITHERED_IMAGE.accum', node='prev')
                           )
    frame = Result(prods.ProcessedImage)
    sky = Result(prods.ProcessedImage, optional=True)
    #
    # Accumulate Frame results
    accum = Result(prods.ProcessedImage, optional=True)

    #@emirdrp.decorators.aggregate
    @emirdrp.decorators.loginfo
    def run(self, rinput):
        partial_result = self.run_single(rinput)
        new_result = self.aggregate_result(partial_result, rinput)
        return new_result

    def aggregate_result(self, partial_result, rinput):

        obresult = rinput.obresult
        # Check if this is our first run
        naccum = getattr(obresult, 'naccum', 1)
        accum = rinput.accum

        frame = partial_result.frame

        if 0 <= naccum <= 1  or accum is None:
            self.logger.debug('round %d initialize accumulator', naccum)
            newaccum = frame
        else:
            self.logger.debug('round %d of accumulation', naccum)
            newaccum = self.aggregate_frames(accum, frame, naccum)

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
        data_arr_sr, regions = resize_arrays(
            [m[0].data for m in data_hdul_s],
            subpixshape,
            offsetsp,
            finalshape,
            fill=1
        )

        if self.intermediate_results:
            self.logger.debug('save resized intermediate img')
            for idx, arr_r in enumerate(data_arr_sr):
                self.save_intermediate_array(arr_r, 'interm1_%03d.fits' % idx)

        hdulist = self.combine(data_arr_sr, data_hdul, finalshape, offsetsp,
                               refpix, use_errors)

        self.save_intermediate_img(hdulist, 'result_initial1.fits')

        compute_cross_offsets = True
        if compute_cross_offsets:

            self.logger.debug("Compute cross-correlation of images")
            # regions = self.compute_regions(finalshape, box=200, corners=True)

            # Regions frm bright objects
            regions = self.compute_regions_from_objs(hdulist[0].data, finalshape, box=20)

            try:

                offsets_xy_c = self.compute_offset_xy_crosscor_regions(
                    data_arr_sr, regions, refine=True, tol=1
                )
    #
                # Combined offsets
                # Offsets in numpy order, swaping
                offsets_xy_t = offset_xy0 - offsets_xy_c
                offsets_fc = offsets_xy_t[:, ::-1]
                offsets_fc_t = numpy.round(offsets_fc).astype('int')
                self.logger.debug('Total offsets: %s', offsets_xy_t)
                self.logger.info('Computing relative offsets from cross-corr')
                finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)
    #
                self.logger.debug("Relative offsetsp (crosscorr) %s", offsetsp)
                self.logger.info('Shape of resized array (crosscorr) is %s', finalshape)

                # Resizing target imgs
                self.logger.debug("Resize to final offsets")
                data_arr_sr, regions = resize_arrays(
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

                hdulist = self.combine(data_arr_sr, data_hdul, finalshape, offsetsp,
                                       refpix, use_errors)

                self.save_intermediate_img(hdulist, 'result_initial2.fits')
            except Exception as error:
                self.logger.warning('Error during cross-correlation, %s', error)

        result = self.create_result(frame=hdulist, sky=sky_result)
        self.logger.info('end of dither recipe')
        return result

    def combine(self, data_arr_sr, data_hdul, finalshape, offsetsp, refpix, use_errors):
        # FIXME: this is mostly duplicated
        # in processing.combine
        baseimg = data_hdul[0]
        has_num_ext = 'NUM' in baseimg
        has_bpm_ext = 'BPM' in baseimg
        baseshape = baseimg[0].shape
        subpixshape = baseshape
        base_header = baseimg[0].header
        result = copy_img(baseimg)

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
        mask_arr_r, _ = resize_arrays(masks, subpixshape, offsetsp, finalshape, fill=1)

        # Position of refpixel in final image
        refpix_final = refpix + offsetsp[0]
        self.logger.info('Position of refpixel in final image %s', refpix_final)

        self.logger.info('Combine target images (final)')
        method = combine.median
        out = method(data_arr_sr, masks=mask_arr_r, dtype='float32')

        self.logger.debug('create result image')
        result[0].data = out[0]
        hdu = result[0]
        self.logger.debug('update result header')
        hdr = hdu.header
        self.set_base_headers(hdr)

        hdr['TSUTC2'] = data_hdul[-1][0].header['TSUTC2']
        # Update obsmode in header

        hdu.header['history'] = "Combined %d images using '%s'" % (
            len(data_hdul),
            method.__name__
        )
        hdu.header['history'] = 'Combination time {}'.format(
            datetime.datetime.utcnow().isoformat()
        )
        # Update NUM-NCOM, sum of individual imagess
        ncom = 0
        for img in data_hdul:
            hdu.header['history'] = "Image {}".format(self.datamodel.get_imgid(img))
            ncom += img[0].header.get('NUM-NCOM', 1)
        hdr['NUM-NCOM'] = ncom
        hdr['UUID'] = str(uuid.uuid1())

        # Update WCS, approximate solution
        hdr['CRPIX1'] += offsetsp[0][0]
        hdr['CRPIX2'] += offsetsp[0][1]
        #
        if use_errors:
            varhdu = fits.ImageHDU(out[1], name='VARIANCE')
            result.append(varhdu)
            num = fits.ImageHDU(out[2].astype('int16'), name='MAP')
            result.append(num)
        return result

    def combine2(self, data, masks, data_hdul, offsetsp, use_errors):
        # FIXME: this is mostly duplicated
        # in processing.combine
        baseimg = data_hdul[0]
        base_header = baseimg[0].header
        result = copy_img(baseimg)

        self.logger.info('Combine target images (final)')
        method = combine.median
        out = method(data, masks=masks, dtype='float32')

        out = quantileclip(data, masks,
                           dtype='float32', out=out, fclip=0.1)

        self.logger.debug('create result image')
        hdu = result[0]
        hdu.data = out[0]
        self.logger.debug('update result header')
        hdr = hdu.header
        self.set_base_headers(hdr)

        hdr['TSUTC2'] = data_hdul[-1][0].header['TSUTC2']
        hdr['UUID'] = str(uuid.uuid1())
        # Update obsmode in header

        hdu.header['history'] = "Combined %d images using '%s'" % (
            len(data_hdul),
            method.__name__
        )
        hdu.header['history'] = 'Combination time {}'.format(
            datetime.datetime.utcnow().isoformat()
        )
        # Update NUM-NCOM, sum of individual images
        ncom = 0
        for img in data_hdul:
            hdu.header['history'] = "Image {}".format(self.datamodel.get_imgid(img))
            ncom += img[0].header.get('NUM-NCOM', 1)
        hdr['NUM-NCOM'] = ncom
        # Update WCS, approximate solution
        hdr['CRPIX1'] += offsetsp[0][0]
        hdr['CRPIX2'] += offsetsp[0][1]

        #
        if use_errors:
            varhdu = fits.ImageHDU(out[1], name='VARIANCE')
            result.append(varhdu)
            num = fits.ImageHDU(out[2].astype('int16'), name='MAP')
            result.append(num)
        return result

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(JoinDitheredImagesRecipe, self).set_base_headers(hdr)
        hdr['IMGOBBL'] = 0
        hdr['OBSMODE'] = 'DITHERED_IMAGE'
        return hdr

    def compute_offset_wcs_imgs(self, imgs, baseshape, subpixshape):

        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2).astype('float')
        offsets_xy = offsets_from_wcs_imgs(imgs, refpix)
        self.logger.debug("offsets_xy %s", offsets_xy)
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:, ::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        self.logger.info('Computing relative offsets')
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)

        return finalshape, offsetsp, refpix, offsets_xy

    def compute_offset_crosscor(self, arrs, region, subpixshape, refine=False):
        offsets_xy = offsets_from_crosscor(arrs, region, refine=refine, order='xy')
        self.logger.debug("offsets_xy cross-corr %s", offsets_xy)
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:, ::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        self.logger.info('Computing relative offsets from cross-corr')
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)

        return finalshape, offsetsp, offsets_xy

    def compute_offset_xy_crosscor_regions(self, arrs, regions, refine=False, tol=0.5):
        offsets_xy = offsets_from_crosscor_regions(
            arrs, regions,
            refine=refine, order='xy', tol=tol
        )
        self.logger.debug("offsets_xy cross-corr %s", offsets_xy)
        # Offsets in numpy order, swaping
        return offsets_xy

    def compute_offset_crosscor_regions(self, arrs, regions, subpixshape, refine=False, tol=0.5):
        offsets_xy = offsets_from_crosscor_regions(
            arrs, regions,
            refine=refine, order='xy', tol=tol
        )
        self.logger.debug("offsets_xy cross-corr %s", offsets_xy)
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:, ::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        self.logger.info('Computing relative offsets from cross-corr')
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)

        return finalshape, offsetsp, offsets_xy


    def compute_shapes_wcs(self, imgs):

        # Better near the center...
        shapes = [img[0].shape for img in imgs]
        ref_pix_xy_0 = (shapes[0][1] // 2, shapes[0][0] // 2)
        #
        ref_coor_xy = reference_pix_from_wcs_imgs(imgs, ref_pix_xy_0)
        # offsets_xy = offsets_from_wcs_imgs(imgs, numpy.asarray([ref_pix_xy_0]))
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

    def compute_sky_simple(self, data_hdul, use_errors=False):

        refimg = data_hdul[0]
        base_header = refimg[0].header
        self.logger.info('combine images with median')
        method = combine.median
        for m in data_hdul:
            m = numpy.median(m[0].data)
            self.logger.debug('median is %f', m)
        sky_data = method([m[0].data for m in data_hdul], dtype='float32')

        hdu = fits.PrimaryHDU(sky_data[0], header=base_header)

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

        if use_errors:
            varhdu = fits.ImageHDU(sky_data[1], name='VARIANCE')
            num = fits.ImageHDU(sky_data[2], name='MAP')
            sky_result = fits.HDUList([hdu, varhdu, num])
        else:
            sky_result = fits.HDUList([hdu])

        return sky_result

    def compute_sky_advanced(self, data_hdul, omasks, base_header, use_errors):
        method = combine.mean

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

    def aggregate2(self, frame1, frame2, naccum):
        # FIXME, this is almost identical to run_single
        frames = [frame1, frame2]
        use_errors = True
        # Initial checks
        fframe = frames[0]
        img = fframe.open()
        result = copy_img(img)
        base_header = result[0].header

        imgs = []
        for f in frames:
            img = f.open()
            imgs.append(img)

        self.logger.info('Computing offsets from WCS information')

        finalshape, partial_shapes, refpix_xy_0, refpix_final_xy = self.compute_shapes_wcs(imgs)

        self.logger.info('Shape of resized array is %s', finalshape)
        self.logger.debug("partial shapes %s", partial_shapes)

        masks = []
        self.logger.debug('Obtains masks')
        for img in imgs:
            if 'NUM' in img:
                self.logger.debug('Using NUM extension as mask')
                mask = numpy.where(img['NUM'].data, 0, 1).astype('int16')
            elif 'BPM' in img:
                self.logger.debug('Using BPM extension as mask')
                mask = numpy.where(img['BPM'].data, 1, 0).astype('int16')
            else:
                self.logger.warning('BPM missing, use zeros instead')
                mask = numpy.zeros_like(img[0].data)
            masks.append(mask)

        # Resizing target frames
        data_arr_r = resize_arrays_alt(
            [img[0].data for img in imgs],
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
        hdu = result[0]
        hdu.data = out[0]
        self.logger.debug('update result header')
        hdr = hdu.header
        self.set_base_headers(hdr)
        hdr['IMGOBBL'] = 0
        hdr['TSUTC2'] = imgs[-1][0].header['TSUTC2']
        # Update obsmode in header
        hdr['OBSMODE'] = 'DITHERED_IMAGE'
        hdu.header['history'] = "Combined %d images using '%s'" % (
            len(imgs),
            method.__name__
        )
        hdu.header['history'] = 'Combination time {}'.format(
            datetime.datetime.utcnow().isoformat()
        )
        # Update NUM-NCOM, sum of individual frames
        ncom = 0
        for img in imgs:
            hdu.header['history'] = "Image {}".format(self.datamodel.get_imgid(img))
            ncom += img[0].header['NUM-NCOM']

        hdr['NUM-NCOM'] = ncom
        # Update WCS, approximate solution
        hdr['CRPIX1'] += (refpix_final_xy[0] - refpix_xy_0[0])
        hdr['CRPIX2'] += (refpix_final_xy[1] - refpix_xy_0[1])

        #
        if use_errors:
            varhdu = fits.ImageHDU(out[1], name='VARIANCE')
            result.append(varhdu)
            num = fits.ImageHDU(out[2].astype('int16'), name='MAP')
            result.append(num)
        return result

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
            region = image_box2d(obj['x'], obj['y'], finalshape, (box, box))
            regions.append(region)
        return regions

    def compute_regions(self, finalshape, box=200, corners=True):
        regions = []
        # A square of 100x100 in the center of the image
        xref_cross = finalshape[1] // 2
        yref_cross = finalshape[0] // 2
        #
        self.logger.debug("Reference position is (x,y) %d  %d", xref_cross + 1, yref_cross + 1)
        self.logger.debug("Reference regions size is %d", 2 * box + 1)
        region = image_box2d(xref_cross, yref_cross, finalshape, (box, box))
        regions.append(region)
        # corners
        if corners:
            xref_c = finalshape[1] // 4
            yref_c = finalshape[0] // 4

            for xi in [xref_c, 3 * xref_c]:
                for yi in [yref_c, 3 * yref_c]:
                    self.logger.debug("Reference position is (x,y) %d  %d", xi + 1, yi + 1)
                    self.logger.debug("Reference regions size is %d", 2 * box + 1)
                    region = image_box2d(xi, yi, finalshape, (box, box))
                    regions.append(region)

        return regions

    def create_object_catalog(self, arr, threshold=1.5, border=0):

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


class FullDitheredImagesRecipe(JoinDitheredImagesRecipe):
    obresult = ObservationResultRequirement(query_opts=ResultOf('frame', node='children'))
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    # extinction = Extinction_Requirement()
    # sources = Catalog_Requirement()
    # offsets = Offsets_Requirement()
    offsets = Requirement(
        prods.CoordinateList2DType,
        'List of pairs of offsets',
        optional=True
    )

    iterations = Parameter(4, 'Iterations of the recipe')
    sky_images = Parameter(
        5, 'Images used to estimate the '
           'background before and after current image')
    sky_images_sep_time = reqs.SkyImageSepTime_Requirement()
    check_photometry_levels = Parameter(
        [0.5, 0.8], 'Levels to check the flux of the objects')
    check_photometry_actions = Parameter(
        ['warn', 'warn', 'default'], 'Actions to take on images')

    frame = Result(prods.ProcessedImage)
    sky = Result(prods.ProcessedImage, optional=True)
    catalog = Result(prods.SourcesCatalog, optional=True)

    def run(self, rinput):
        partial_result = self.run_single(rinput)
        return partial_result

    def run_single(self, rinput):

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
        data_arr_sr, regions = resize_arrays(
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
        mask_arr_r, _ = resize_arrays(masks, subpixshape, offsetsp, finalshape, fill=1)

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
                finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)
                #
                self.logger.debug("Relative offsetsp (crosscorr) %s", offsetsp)
                self.logger.info('Shape of resized array (crosscorr) is %s', finalshape)

                # Resizing target imgs
                self.logger.debug("Resize to final offsets")
                data_arr_sr, regions = resize_arrays(
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
                mask_arr_r, _ = resize_arrays(masks, subpixshape, offsetsp, finalshape, fill=1)

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
                sky, _, num = median(skydata, skymasks, scales=skyscales)
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
                    fixpix2(sky, binmask, out=sky, iterations=1)

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

    def compute_superflat(self, data_arr_r, objmask, regions, channels):
        # superflat

        mask = [objmask[r] for r in regions]
        scales = [numpy.median(d) for d in data_arr_r]
        self.logger.debug('flat scaling %s', scales)
        sf_data, _sf_var, sf_num = flatcombine(data_arr_r, masks=mask, scales=scales)

        for channel in channels:
            mask = (sf_num[channel] == 0)
            if numpy.any(mask):
                fixpix2(sf_data[channel], mask, out=sf_data[channel])

        # Normalize, flat has mean = 1
        sf_data /= sf_data.mean()
        return sf_data

    def compute_sky_advanced(self, data_hdul, omasks, base_header, use_errors):
        method = combine.mean

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
