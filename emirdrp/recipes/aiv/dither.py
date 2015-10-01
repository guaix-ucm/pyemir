#
# Copyright 2014-2015 Universidad Complutense de Madrid
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

'''AIV Recipes for EMIR'''

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

from emirdrp.core import offsets_from_wcs
from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType


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


def combine_frames(rframes):
    # frameslll = [frame.open() for frame in rframes]
    frameslll = rframes
    data = [i[0].data for i in frameslll]
    out = combine.mean(data, dtype='float32')
    return out


class DitheredImageARecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    frame = Product(DataFrameType)

    def run(self, rinput):

        _logger.info('Computing offsets from WCS information')
        baseshape = (2048, 2048)
        refpix = numpy.divide(numpy.array([baseshape],
                                          dtype='int'), 2).astype('float')
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
        out = combine_frames(rframes)
        hdu = fits.PrimaryHDU(out[0])

        _logger.debug('update result header')
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdr['IMGOBBL'] = 0

        hdulist = fits.HDUList([hdu])

        result = self.create_result(frame=hdulist)

        return result


from numina.core import ObservationResult


class DitheredImageRecipeInputBuilder(object):
    '''Class to build DitheredImageRecipe inputs from the Observation Results
   
    RecipeInputBuilder which fetches the pre-reduced images that will be combined
   
    '''

    def __init__(self, dal):
        self.dal = dal
   
    def buildRecipeInput(self, obsres):
       
        stareImages = []

        stareImagesIds = obsres['stareImagesIds']._v 
        for subresId in stareImagesIds:
            subres = self.dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['frame'])
        
        newOR = ObservationResult()
        newOR.frames = stareImages
        obsres['obresult'] = newOR
        print('Adding RI parameters ', obsres)
        newRI = DitheredImageARecipeInput(**obsres)

        return newRI


DitheredImageARecipeInput = DitheredImageARecipe.RecipeInput
DitheredImageARecipe.InputBuilder = DitheredImageRecipeInputBuilder

