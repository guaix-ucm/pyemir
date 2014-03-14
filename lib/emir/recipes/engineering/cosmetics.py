#
# Copyright 2010-2014 Universidad Complutense de Madrid
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
from astropy.io import fits

from numina.core import BaseRecipe, Parameter, DataFrame
from numina.core import RecipeError,RecipeRequirements
from numina.core import Product, define_requirements, define_result
from numina.core import FrameDataProduct
from numina.array.cosmetics import cosmetics, PIXEL_DEAD, PIXEL_HOT
from numina.core.requirements import ObservationResultRequirement
from numina.core.requirements import InstrumentConfigurationRequirement

from emir.core import RecipeResult

_logger = logging.getLogger('numina.recipes.emir')

class CosmeticsRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    insconf = InstrumentConfigurationRequirement()
    lowercut = Parameter(4.0, 'Values bellow this sigma level are flagged as dead pixels')
    uppercut = Parameter(4.0, 'Values above this sigma level are flagged as hot pixels')
    maxiter = Parameter(30, 'Maximum number of iterations')
    
class CosmeticsRecipeResult(RecipeResult):
    ratio = Product(FrameDataProduct)
    mask = Product(FrameDataProduct)
    
@define_requirements(CosmeticsRecipeRequirements)
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

    def run(self, rinput):
        
        # FIXME:
        # We need 2 flats
        # Of different exposure times
        #
        # And their calibrations
        #
        if len(rinput.obresult.frames) < 2:
            raise RecipeError('The recipe requires 2 flat frames')
        
        with fits.open(rinput.obresult.frames[0]) as hdul:
            f1 = hdul[0].data.copy()

        with fits.open(rinput.obresult.frames[1]) as hdul:
            f2 = hdul[0].data.copy()

        maxiter = rinput.maxiter
        lowercut = rinput.lowercut
        uppercut = rinput.uppercut
        
        ninvalid = 0
        mask = None
        
        if mask:
            m = fits.getdata(mask)
            ninvalid = numpy.count_nonzero(m)
        else:
            m = numpy.zeros_like(f1, dtype='int')
    
        for niter in range(1, maxiter + 1):
            _logger.info('iter %d', niter)
            ratio, m, sigma = cosmetics(f1, f2, m, lowercut=lowercut, uppercut=uppercut)
            # FIXME
            # These are intermediate results that 
            # can be removed later
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fits.writeto('numina-cosmetics-i%02d.fits' % niter, ratio, clobber=True)
                fits.writeto('numina-mask-i%02d.fits' % niter, m, clobber=True)
                fits.writeto('numina-sigma-i%02d.fits' % niter, m * 0.0 + sigma, clobber=True)
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
            fits.writeto('numina-cosmetics.fits', ratio, clobber=True)
            fits.writeto('numina-mask.fits', m, clobber=True)
            fits.writeto('numina-sigma.fits', sigma * numpy.ones_like(m), clobber=True)
            
        hdu = fits.PrimaryHDU(ratio)
        hdr = hdu.header
        hdr['NUMXVER'] = ( __version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        ratiohdl = fits.HDUList([hdu])    
        
        maskhdu = fits.PrimaryHDU(m)
        hdr = maskhdu.header
        hdr['NUMXVER'] = ( __version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        maskhdl = fits.HDUList([maskhdu])

        res = CosmeticsRecipeResult(ratio=ratiohdl, mask=maskhdl)
        return res

