#
# Copyright 2010-2012 Universidad Complutense de Madrid
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

'''Recipe for finding cosmetic defects.'''

from __future__ import division

import logging
import warnings

import numpy
import pyfits

from numina.core import BaseRecipe, Parameter, DataFrame
from numina.core import RecipeError,RecipeInput, RecipeResult
from numina.core import Product, define_input, define_result
from numina.core import FrameDataProduct
from numina.array.cosmetics import cosmetics, PIXEL_DEAD, PIXEL_HOT

_logger = logging.getLogger('numina.recipes.emir')

class CosmeticsRecipeInput(RecipeInput):
    lowercut = Parameter(4.0, 'Values bellow this sigma level are flagged as dead pixels')
    uppercut = Parameter(4.0, 'Values above this sigma level are flagged as hot pixels')
    maxiter = Parameter(30, 'Maximum number of iterations')
    
class CosmeticsRecipeResult(RecipeResult):
    ratio = Product(FrameDataProduct)
    mask = Product(FrameDataProduct)
    
@define_input(CosmeticsRecipeInput)
@define_result(CosmeticsRecipeResult)
class CosmeticsRecipe(BaseRecipe):
    '''Detector Cosmetics.
    
    Recipe to find and tag bad pixels in the detector.
    '''

    def __init__(self):
        super(CosmeticsRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):

        resets = []
        flats = []

        for frame in obresult.frames:
            if frame.itype == 'RESET':
                resets.append(frame.label)
                _logger.debug('%s is RESET', frame.label)
            elif frame.itype == 'FLAT':
                flats.append(frame.label)
                _logger.debug('%s is FLAT', frame.label)
            else:
                raise RecipeError('frame is neither a FLAT nor a RESET')


        # we need 2 flats and 1 reset
        if len(flats) < 2:
            raise ValueError('The recipe requires 2 flat frames')
        
        if len(resets) < 1:
            raise ValueError('The recipe requires 1 reset frame')
        
        reset = pyfits.getdata(resets[-1])
        f1 = pyfits.getdata(flats[0]) - reset
        f2 = pyfits.getdata(flats[1]) - reset
        
        maxiter = self.parameters['maxiter']
        lowercut = self.parameters['lowercut']
        uppercut = self.parameters['uppercut']
        
        ninvalid = 0
        mask = None
        
        if mask:
            m = pyfits.getdata(mask)
            ninvalid = numpy.count_nonzero(m)
        else:
            m = numpy.zeros_like(reset, dtype='int')
    
        for niter in range(1, maxiter + 1):
            _logger.info('iter %d', niter)
            ratio, m, sigma = cosmetics(f1, f2, m, lowercut=lowercut, uppercut=uppercut)
            # FIXME
            # These are intermediate results that 
            # can be removed later
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pyfits.writeto('numina-cosmetics-i%02d.fits' % niter, ratio, clobber=True)
                pyfits.writeto('numina-mask-i%02d.fits' % niter, m, clobber=True)
                pyfits.writeto('numina-sigma-i%02d.fits' % niter, m * 0.0 + sigma, clobber=True)
            _logger.info('iter %d, invalid points in input mask: %d', niter, ninvalid)
            _logger.info('iter %d, estimated sigma is %f', niter, sigma)
            n_ninvalid = numpy.count_nonzero(m)
            
            # Probably there is something wrong here
            # too much defective pixels
            if ninvalid / m.size >= 0.10:
                # This should set a flag in the output
                msg = 'defective pixels are greater than 10%'
                _logger.warning(msg)
                            
            if n_ninvalid == ninvalid:
                _logger.info('convergence reached after %d iterations', niter)
                break
            _logger.info('new invalid points: %d', n_ninvalid - ninvalid)            
            ninvalid = n_ninvalid
        else:
            # This should set a flag in the output
            msg = 'convergence not reached after %d iterations' % maxiter
            _logger.warning(msg)
    
        _logger.info('number of dead pixels %d', numpy.count_nonzero(m == PIXEL_DEAD))
        _logger.info('number of hot pixels %d', numpy.count_nonzero(m == PIXEL_HOT))
    
        # FIXME
        # These are intermediate results that 
        # caN be removed later
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pyfits.writeto('numina-cosmetics.fits', ratio, clobber=True)
            pyfits.writeto('numina-mask.fits', m, clobber=True)
            pyfits.writeto('numina-sigma.fits', sigma * numpy.ones_like(m), clobber=True)
            
        ratiohdu = pyfits.PrimaryHDU(ratio)
        hdr = ratiohdu.header
        hdr.update('FILENAME', 'ratio.fits')
        hdr = self.update_header(hdr)        
        # hdr.update('IMGTYP', 'TARGET', 'Image type')
        # hdr.update('NUMTYP', 'TARGET', 'Data product type')
        ratiohdl = pyfits.HDUList([ratiohdu])    
        
        maskhdu = pyfits.PrimaryHDU(m)
        hdr = maskhdu.header
        hdr.update('filename', 'mask.fits')
        hdr = self.update_header(hdr)
        maskhdl = pyfits.HDUList([maskhdu])
        

        res = CosmeticsRecipeResult(ratio=DataFrame(ratiohdl), 
                                        mask=DataFrame(maskhdl))        
        return res
