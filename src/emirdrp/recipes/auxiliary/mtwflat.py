#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""Twilight Flat Recipe for a list of frames in different filters"""


import uuid
import datetime

import numpy
import astropy.io.fits as fits
from numina.array.combine import median
from numina.array.robustfit import fit_theil_sen
from numina.core import Result
from numina.frame.utils import copy_img
import numina.types.datatype as dt
from numina.processing.combine import basic_processing_with_combination_frames

from emirdrp.processing.info import gather_info_frames
from emirdrp.core.recipe import EmirRecipe
import emirdrp.products as prods
import emirdrp.requirements as reqs


class MultiTwilightFlatRecipe(EmirRecipe):
    """Create a list of twilight flats"""
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    twflatframes = Result(dt.ListOfType(prods.MasterIntensityFlat))

    def run(self, rinput):

        results = []
        self.logger.info('starting multiflat flat reduction')

        # Uncomment this line
        # to revert to non-ramp
        # flow = self.init_filters(rinput)
        saturation = 45000.0

        iinfo = gather_info_frames(rinput.obresult.frames)
        image_groups = {}
        self.logger.info('group images by filter')
        for idx, info in enumerate(iinfo):
            filt = info['filter']
            if filt not in image_groups:
                self.logger.debug('new filter %s', filt)
                image_groups[filt] = []
            img = rinput.obresult.frames[idx]
            self.logger.debug('image %s in group %s', img, filt)
            image_groups[filt].append(img)

        for filt, frames in image_groups.items():
            self.logger.info('processing filter %s', filt)

            # Uncomment this line and comment the following
            # to revert to non-ramp
            # res = self.run_per_filter(frames, flow)
            try:
                res = self.run_per_filter_ramp(frames, saturation=saturation)
                results.append(res)
            except ValueError:
                self.logger.info('filter %s cannot be processed', filt)

        self.logger.info('end multiflat flat reduction')
        result = self.create_result(twflatframes=results)

        return result

    def run_per_filter(self, frames, flow):

        errors = True
        self.logger.debug('using errors: %s', errors)
        hdulist = basic_processing_with_combination_frames(frames, flow,
                                                    method=median,
                                                    errors=errors)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        mm = hdulist[0].data.mean()
        self.logger.info('mean value of flat is: %f', mm)
        hdr['CCDMEAN'] = mm

        self.logger.debug('normalize image')
        hdulist[0].data /= mm
        if errors:
            self.logger.debug('normalize VAR extension')
            hdulist['variance'].data /= (mm * mm)

        return hdulist

    def run_per_filter_ramp(self, frames, saturation, errors=False):
        imgs = [frame.open() for frame in frames]
        return self.run_img_per_filter_ramp(imgs, saturation, errors)

    def run_img_per_filter_ramp(self, imgs, saturation, errors=False):

        nimages = len(imgs)
        if nimages == 0:
            raise ValueError('len(images) == 0')

        median_frames = numpy.empty((nimages,))
        exptime_frames = []
        utc_frames = []
        # generate 3D cube
        bshape = self.datamodel.shape
        flat_frames = numpy.empty((bshape[0], bshape[1], nimages))
        for idx, image in enumerate(imgs):
            flat_frames[:, :, idx] = image['primary'].data
            exptime_frames.append(image[0].header['EXPTIME'])
            median_frames[idx] = numpy.median(image['primary'].data)
            utc_frames.append(image[0].header['UTC'])

            self.logger.debug(
                'image %d exptime %f median %f UTC %s',
                idx, exptime_frames[idx],
                median_frames[idx],
                utc_frames[-1]
            )

        # filter saturated images
        good_images = median_frames < saturation
        ngood_images = good_images.sum()
        slope_scaled_var = None
        slope_scaled_num = None
        if ngood_images == 0:
            self.logger.warning('We have only %d good images', ngood_images)
            raise ValueError('No images under saturation')
        elif ngood_images < 2:
            self.logger.warning('We have only %d good images', ngood_images)
            # Reference image
            ref_image = imgs[0]
            slope_scaled = numpy.ones(bshape) * exptime_frames[0]

            if errors:
                slope_scaled_var = numpy.zeros_like(slope_scaled)
                slope_scaled_num = numpy.zeros_like(slope_scaled, dtype='int16') + ngood_images
        else:
            nsaturated = nimages - good_images.sum()
            if nsaturated > 0:
                self.logger.debug('we have %d images with median value over saturation (%f)',
                                  nsaturated , saturation)

            m = flat_frames[:,:, good_images]
            # Reshape array to obtain a 2D array
            m_r = m.reshape((bshape[0] * bshape[1], ngood_images))
            self.logger.debug('fitting slopes with Theil-Sen')
            # self.logger.debug('fitting slopes with mean-squares')
            # ll = nppol.polyfit(median_frames[good_images], m_r.T, deg=1)
            ll = self.filter_nsigma(median_frames[good_images], m_r.T)

            slope = ll[1].reshape(bshape)
            base = ll[0].reshape(bshape)

            # First good frame
            index_of_first_good = numpy.nonzero(good_images)[0][0]
            slope_scaled = slope * exptime_frames[index_of_first_good]
            if errors:
                slope_scaled_var = numpy.zeros_like(slope_scaled)
                slope_scaled_num = numpy.zeros_like(slope_scaled, dtype='int16') + ngood_images

        cdata = []
        for idx, img in enumerate(imgs):
            if good_images[idx]:
                cdata.append(img)

        result = self.compose_result(cdata, slope_scaled, errors, slope_scaled_var, slope_scaled_num)

        return result

    def filter_nsigma(self, median_val, image_val, nsigma=10.0, nloop=1):
        # Initial estimation
        ll = fit_theil_sen(median_val, image_val)
        ni = 0
        self.logger.debug('initial estimation')
        while ni < nloop:
            # Prediction
            self.logger.debug('loop %d', ni + 1)
            base, slope = ll
            image_val_pred = base + median_val[:,numpy.newaxis] * slope
            image_diff = image_val - image_val_pred
            # Compute MAD
            mad = compute_mad(image_diff)
            sigma_robust = nsigma * 1.4826 * mad
            self.logger.debug('compute robust std deviation')
            self.logger.debug(
                'min %7.1f max %7.1f mean %7.1f',
                sigma_robust.min(),
                sigma_robust.max(),
                sigma_robust.mean()
            )
            # Check values over sigma
            mask_over = numpy.abs(image_diff) >= sigma_robust[:, numpy.newaxis]
            self.logger.debug('values over sigma: %d', mask_over.sum())
            # Insert expected values in image
            # instead of masking
            image_val[mask_over] = image_val_pred[mask_over]
            #
            self.logger.debug('Theil-Sen fit')
            ll = fit_theil_sen(median_val, image_val)
            ni += 1

        return ll

    def compose_result(self, imgs, slope_scaled, errors=False, slope_scaled_var=None, slope_scaled_num=None):

        self.logger.debug('update result header')
        cnum = len(imgs)
        method_name = 'Theil-Sen'
        result = copy_img(imgs[0])
        base_header = result[0].header
        cdata = imgs

        hdu = result[0]
        hdu.data = slope_scaled
        self.set_base_headers(hdu.header)

        hdu.header['history'] = "Combined %d images using '%s'" % (cnum, method_name)
        hdu.header['history'] = 'Combination time {}'.format(datetime.datetime.utcnow().isoformat())

        for img in cdata:
            hdu.header['history'] = "Image {}".format(self.datamodel.get_imgid(img))

        prevnum = base_header.get('NUM-NCOM', 1)
        hdu.header['NUM-NCOM'] = prevnum * cnum
        hdu.header['UUID'] = str(uuid.uuid1())
        # Headers of last image
        hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
        # TODO: use BPM to compute mean
        mm = hdu.data.mean()
        self.logger.info('mean value of flat is: %f', mm)
        hdu.header['CCDMEAN'] = mm
        self.logger.debug('normalize image')
        hdu.data /= mm

        if errors:
            varhdu = fits.ImageHDU(slope_scaled_var, name='VARIANCE')
            result.append(varhdu)
            num = fits.ImageHDU(slope_scaled_num, name='MAP')
            result.append(num)
            self.logger.debug('normalize VAR extension')
            varhdu.data /= (mm * mm)
        return result


def compute_mad(x):
    m1 = numpy.median(x, axis=1)
    m2 = numpy.abs(x - m1[:, numpy.newaxis])
    mad = numpy.median(m2, axis=1)
    return mad
