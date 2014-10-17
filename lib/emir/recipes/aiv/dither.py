#
# Copyright 2014 Universidad Complutense de Madrid
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

from __future__ import division

import logging

import numpy
from astropy.io import fits


from numina import __version__
from numina.core import BaseRecipe, RecipeRequirements
from numina.core import Product
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement
from numina.array import combine
from numina.array import combine_shape
from numina.array import subarray_match
from numina.frame import resize_hdu

from emir.core import offsets_from_wcs
from emir.core import RecipeResult
from emir.core import EMIR_BIAS_MODES
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import DataFrameType, MasterIntensityFlat
from emir.dataproducts import CoordinateList2DType
from emir.dataproducts import ArrayType
from emir.core import gather_info

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"
          

def resize_hdul(hdul, newshape, region, extensions=None, window=None,
                scale=1, fill=0.0, clobber=True, conserve=True):
           
    if extensions is None:
        extensions = [0]
    
    nhdul = [None] * len(hdul)
    for ext, hdu in enumerate(hdul):
        if ext in extensions:
            nhdul[ext] = resize_hdu(hdu, newshape, region, fill=fill, window=window, scale=scale, conserve=conserve)
        else:
            nhdul[ext] = hdu
    return fits.HDUList(nhdul)
    
       
def resize(frames, shape, offsetsp, finalshape, window=None, scale=1, step=0):
    _logger.info('Resizing frames and masks')
    rframes = []
    for frame, rel_offset in zip(frames, offsetsp):
        region, _ = subarray_match(finalshape, rel_offset, shape)
        print frame, region
        rframe = resize_hdul(frame.open(), finalshape, region)
        rframes.append(rframe)
                
    return rframes
    
def combine_frames(rframes):
    #frameslll = [frame.open() for frame in rframes]
    frameslll = rframes
    data = [i[0].data for i in frameslll]
    out = combine.mean(data, dtype='float32')
    return out
            
class DitheredImageARecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()

class DitheredImageARecipeResult(RecipeResult):
    frame = Product(DataFrameType)

@define_requirements(DitheredImageARecipeRequirements)
@define_result(DitheredImageARecipeResult)
class DitheredImageARecipe(BaseRecipe):

    def __init__(self):
        super(DitheredImageARecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
                
        _logger.info('Computing offsets from WCS information')
        baseshape = (2048, 2048)
        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2).astype('float')
        offsets_xy = offsets_from_wcs(rinput.obresult.frames, refpix)
        print offsets_xy
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:,::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        _logger.info('Computing relative offsets')
        subpixshape = (2048, 2048)
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)
        print offsetsp
        _logger.info('Shape of resized array is %s', finalshape)
                
        # Resizing target frames 
        rframes = resize(rinput.obresult.frames, subpixshape, offsetsp, finalshape)
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

