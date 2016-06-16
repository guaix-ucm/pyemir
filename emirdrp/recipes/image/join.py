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

import numpy
from astropy.io import fits
from numina import __version__
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
from numina.array import combine
from numina.array import combine_shape
from numina.array import subarray_match
from numina.frame import resize_hdu
from numina.core import ObservationResult

from emirdrp.core import offsets_from_wcs
from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.ext.gtc import RUN_IN_GTC

_logger = logging.getLogger('numina.recipes.emir')


def resize_hdul(hdul, newshape, region, extensions=None, window=None,
                scale=1, fill=0.0, clobber=True, conserve=True):

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
    _logger.info('Resizing frames and masks')
    rframes = []
    for frame, rel_offset in zip(frames, offsetsp):
        region, _ = subarray_match(finalshape, rel_offset, shape)
        rframe = resize_hdul(frame.open(), finalshape, region)
        rframes.append(rframe)

    return rframes


def resize_arrays(arrays, shape, offsetsp, finalshape, window=None, scale=1, conserve=True):

    _logger.info('Resizing arrays')
    rarrays = []
    for array, rel_offset in zip(arrays, offsetsp):
        region, _ = subarray_match(finalshape, rel_offset, shape)
        newdata = resize_array(array, finalshape, region, window=window,
                               fill=0.0, scale=scale, conserve=conserve)
        rarrays.append(newdata)

    return rarrays


def resize_array(data, finalshape, region, window=None,
                     scale=1, fill=0.0, conserve=True):
    from numina.array import rebin_scale
    if window:
        data = data[window]

    if scale == 1:
        finaldata = data
    else:
        finaldata = rebin_scale(data, scale)

    newdata = numpy.empty(finalshape, dtype='float')
    newdata.fill(fill)
    newdata[region] = finaldata
    # Conserve the total sum of the original data
    if conserve:
        newdata[region] /= scale ** 2
    return newdata


def combine_frames(rframes):
    # frameslll = [frame.open() for frame in rframes]
    frameslll = rframes
    data = [i[0].data for i in frameslll]
    out = combine.mean(data, dtype='float32')
    return out


class JoinDitheredImagesRecipe(EmirRecipe):
    """Combine single exposures obtained in dithered mode"""

    obresult = ObservationResultRequirement()
    frame = Product(DataFrameType)

    @classmethod
    def build_recipe_input(cls, obsres, dal, pipeline='default'):
        if RUN_IN_GTC:
            _logger.debug('Using GTC version of build_recipe_input in DitheredImages')
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

        fframe = rinput.obresult.frames[0]
        img = fframe.open()
        template_header = img[0].header

        _logger.debug('Data types of input images')
        data_arrays = []
        for f in rinput.obresult.frames:
            img = f.open()
            _logger.debug('datatype is %s', img[0].data.dtype)
            data_arrays.append(img[0].data)

        _logger.info('Computing offsets from WCS information')
        baseshape = (2048, 2048)
        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2).astype('float')
        offsets_xy = offsets_from_wcs(rinput.obresult.frames, refpix)
        _logger.debug("offsets_xy %s", offsets_xy)
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:, ::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        _logger.info('Computing relative offsets')
        subpixshape = (2048, 2048)
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)
        _logger.debug("offsetsp %s", offsetsp)
        _logger.info('Shape of resized array is %s', finalshape)

        # Resizing target frames
        rframes = resize(rinput.obresult.frames, subpixshape,
                         offsetsp, finalshape)
        r_arrays = resize_arrays(data_arrays, subpixshape, offsetsp, finalshape)
        for f in r_arrays:
            _logger.debug('datatype is %s', f.dtype)

        out = combine.mean(r_arrays, dtype='float32')

        hdu = fits.PrimaryHDU(out[0], header=template_header)

        _logger.debug('update result header')
        hdr = hdu.header
        hdr['IMGOBBL'] = 0

        hdulist = fits.HDUList([hdu])

        result = self.create_result(frame=hdulist)

        return result
