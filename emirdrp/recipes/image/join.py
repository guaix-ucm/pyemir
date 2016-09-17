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

"""Image recipes for EMIR for EMIR"""

from __future__ import division, print_function

import logging

import uuid
import numpy
from astropy.io import fits
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
from numina.array import combine
from numina.array import combine_shape
from numina.array import resize_array
from numina.core import ObservationResult
from numina.flow.processing import SkyCorrector

from emirdrp.processing.wcs import offsets_from_wcs
from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.ext.gtc import RUN_IN_GTC
from emirdrp.processing.combine import resize
from emirdrp.processing.combine import segmentation_combined
from emirdrp.processing.datamodel import EmirDataModel


_logger = logging.getLogger('numina.recipes.emir')


def resize_hdul(hdul, newshape, region, extensions=None, window=None,
                scale=1, fill=0.0, conserve=True):
    from numina.frame import resize_hdu
    if extensions is None:
        extensions = [0]

    nhdul = [None] * len(hdul)
    for ext, hdu in enumerate(hdul):
        if ext in extensions:
            nhdul[ext] = resize_hdu(hdu, newshape,
                                    region, fill=fill,
                                    window=window,
                                    scale=scale,
                                    conserve=conserve)
        else:
            nhdul[ext] = hdu
    return fits.HDUList(nhdul)



def resize(frames, shape, offsetsp, finalshape, window=None, scale=1, step=0):
    from numina.array import subarray_match
    _logger.info('Resizing frames and masks')
    rframes = []
    regions = []
    for frame, rel_offset in zip(frames, offsetsp):
        region, _ = subarray_match(finalshape, rel_offset, shape)
        rframe = resize_hdul(frame.open(), finalshape, region)
        rframes.append(rframe)
        regions.append(region)

    return rframes, regions


def resize_arrays(arrays, shape, offsetsp, finalshape, window=None, scale=1, conserve=True, fill=0.0):
    from numina.array import subarray_match
    _logger.info('Resizing arrays')
    rarrays = []
    regions = []
    for array, rel_offset in zip(arrays, offsetsp):
        region, _ = subarray_match(finalshape, rel_offset, shape)
        newdata = resize_array(array, finalshape, region, window=window,
                               fill=fill, scale=scale, conserve=conserve)
        rarrays.append(newdata)
        regions.append(region)
    return rarrays, regions


def resize_array(data, finalshape, region, window=None,
                 scale=1, fill=0.0, conserve=True):
    from numina.array import rebin_scale
    if window is not None:
        data = data[window]

    if scale == 1:
        finaldata = data
    else:
        finaldata = rebin_scale(data, scale)

    newdata = numpy.empty(finalshape, dtype=data.dtype)
    newdata.fill(fill)
    newdata[region] = finaldata
    # Conserve the total sum of the original data
    if conserve:
        newdata[region] /= scale ** 2
    return newdata


class JoinDitheredImagesRecipe(EmirRecipe):
    """Combine single exposures obtained in dithered mode"""

    obresult = ObservationResultRequirement()
    frame = Product(DataFrameType)
    sky = Product(DataFrameType, optional=True)

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
        print('STARE IMAGES IDS: ', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['frame'])

        newOR = ObservationResult()
        newOR.frames = stareImages
        # obsres['obresult'] = newOR
        # print('Adding RI parameters ', obsres)
        # newRI = DitheredImageARecipeInput(**obsres)
        newRI = cls.create_input(obresult=newOR)
        return newRI

    def run(self, rinput):

        use_errors = True
        fframe = rinput.obresult.frames[0]
        img = fframe.open()
        base_header = img[0].header
        baseshape = img[0].data.shape

        if 'NUM-SK' not in base_header:
            compute_sky = True
        else:
            compute_sky = False

        data_hdul = []
        for f in rinput.obresult.frames:
            img = f.open()
            data_hdul.append(img)

        self.logger.info('Computing offsets from WCS information %d', len(data_hdul))

        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2).astype('float')
        offsets_xy = offsets_from_wcs(rinput.obresult.frames, refpix)
        self.logger.debug("offsets_xy %s", offsets_xy)
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:, ::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        self.logger.info('Computing relative offsets')
        subpixshape = baseshape
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)
        self.logger.debug("Relative offsetsp %s", offsetsp)
        self.logger.info('Shape of resized array is %s', finalshape)

        # Resizing target frames

        data_arr_r, regions = resize_arrays([m[0].data for m in data_hdul], subpixshape, offsetsp, finalshape, fill=1)
        self.logger.warning('BPM missing, use zeros instead')
        false_mask = numpy.zeros(baseshape, dtype='int16')
        self.logger.debug('resize bad pixel masks')
        mask_arr_r, _ = resize_arrays([false_mask for _ in data_arr_r], subpixshape, offsetsp, finalshape, fill=1)
        skyid = None
        if compute_sky:
            method = combine.mean
            bpm = None
            self.logger.info("stacking %d images, with offsets using '%s'", len(data_arr_r), method.func_name)
            data1 = method(data_arr_r, masks=mask_arr_r, dtype='float32')

            segmap = segmentation_combined(data1[0])
            # submasks
            if bpm is None:
                omasks = [(segmap[region] > 0) for region in regions]
            else:
                omasks = [((segmap[region] > 0) & bpm) for region in regions]

            self.logger.info("stacking %d images, with objects mask using '%s'", len(data_arr_r), method.func_name)

            sky_data = method([m[0].data for m in data_hdul], masks=omasks, dtype='float32')
            hdu = fits.PrimaryHDU(sky_data[0], header=base_header)
            points_no_data = (sky_data[2] == 0).sum()

            self.logger.debug('update result header')
            skyid= uuid.uuid1().hex
            hdu.header['EMIRUUID'] = skyid
            hdu.header['history'] = "Combined %d images using '%s'" % (len(data_arr_r), method.func_name)
            self.logger.info("missing points, total: %d, fraction: %3.1f", points_no_data, points_no_data / sky_data[2].size)

            if use_errors:
                varhdu = fits.ImageHDU(sky_data[1], name='VARIANCE')
                num = fits.ImageHDU(sky_data[2], name='MAP')
                sky_result = fits.HDUList([hdu, varhdu, num])
            else:
                sky_result = fits.HDUList([hdu])

            self.logger.info('sky correction in individual images')
            data_arr_s = [m[0].data - sky_data[0] for m in data_hdul]
            data_arr_sr, _ = resize_arrays([f for f in data_arr_s], subpixshape, offsetsp, finalshape, fill=1)
        else:
            sky_result = None
            data_arr_sr = data_arr_r

        # Position of refpixel in final image
        refpix_final = refpix + offsetsp[0]
        self.logger.info('Position of refpixel in final image %s', refpix_final)
        out = combine.mean(data_arr_sr, masks=mask_arr_r, dtype='float32')

        hdu = fits.PrimaryHDU(out[0], header=base_header)
        self.logger.debug('update result header')
        hdr = hdu.header
        self.set_base_headers(hdr)
        hdr['IMGOBBL'] = 0
        hdr['TSUTC2'] = data_hdul[-1][0].header['TSUTC2']
        hdr['OBSMODE'] = 'DITHERED_IMAGE'
        hdu.header['history'] = "Combined %d images using '%s'" % (len(data_arr_sr), 'mean')
        # Approximate solution
        hdr['CRPIX1'] += offsetsp[0][0]
        hdr['CRPIX2'] += offsetsp[0][1]

        if compute_sky:
            hdr['NUM-SK'] = skyid
            hdr['history'] = 'Sky subtraction with {}'.format(skyid)
        #
        if use_errors:
            varhdu = fits.ImageHDU(out[1], name='VARIANCE')
            num = fits.ImageHDU(out[2], name='MAP')
            hdulist = fits.HDUList([hdu, varhdu, num])
        else:
            hdulist = fits.HDUList([hdu])

        result = self.create_result(frame=hdulist, sky=sky_result)

        return result
