#
# Copyright 2010-2019 Universidad Complutense de Madrid
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

"""Recipe for finding cosmetic defects."""

from __future__ import division

import logging
import warnings

import numpy
from astropy.io import fits

from numina import __version__
from numina.core import Parameter
from numina.exceptions import RecipeError
from numina.core import Result
from numina.array.cosmetics import cosmetics, PIXEL_DEAD, PIXEL_HOT
from numina.core.requirements import ObservationResultRequirement
import numina.util.flow
import numina.util.node
import numina.processing as proc

from emirdrp.core.recipe import EmirRecipe
import emirdrp.products as prods
import emirdrp.requirements as reqs


_logger = logging.getLogger(__name__)


def gather_info(hdulist):
    n_ext = len(hdulist)

    # READMODE is STRING
    readmode = hdulist[0].header.get('READMODE', 'undefined')
    bunit = hdulist[0].header.get('BUNIT', 'undefined')
    texp = hdulist[0].header.get('EXPTIME')
    adu_s = True
    if bunit:
        if bunit.lower() == 'adu':
            adu_s = False
        elif bunit.lower() == 'adu/s':
            adu_s = True
        else:
            _logger.warning('Unrecognized value for BUNIT %s', bunit)

    return {'n_ext': n_ext,
            'readmode': readmode,
            'texp': texp,
            'adu_s': adu_s}


class CosmeticsRecipe(EmirRecipe):

    """Detector Cosmetics.

    Recipe to find and tag bad pixels in the detector.
    """

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    lowercut = Parameter(
        4.0, 'Values below this sigma level are flagged as dead pixels')
    uppercut = Parameter(
        4.0, 'Values above this sigma level are flagged as hot pixels')
    maxiter = Parameter(30, 'Maximum number of iterations')

    ratioframe = Result(prods.ProcessedImage)
    maskframe = Result(prods.MasterBadPixelMask)

    def run(self, rinput):

        # FIXME:
        # We need 2 flats
        # Of different exposure times
        #
        # And their calibrations
        #

        if len(rinput.obresult.frames) < 2:
            raise RecipeError('The recipe requires 2 flat frames')

        iinfo = []
        for frame in rinput.obresult.frames:
            with frame.open() as hdulist:
                iinfo.append(gather_info(hdulist))

        # Loading calibrations
        with rinput.master_bias.open() as hdul:
            readmode = hdul[0].header.get('READMODE', 'undefined')
            if readmode.lower() in ['simple', 'bias']:
                self.logger.debug('loading bias')
                mbias = hdul[0].data
                bias_corrector = proc.BiasCorrector(mbias)
            else:
                self.logger.debug('ignoring bias')
                bias_corrector = numina.util.node.IdNode()

        with rinput.master_dark.open() as mdark_hdul:
            self.logger.debug('loading dark')
            mdark = mdark_hdul[0].data
            dark_corrector = proc.DarkCorrector(mdark)

        flow = numina.util.flow.SerialFlow([bias_corrector, dark_corrector])

        self.logger.info('processing flat #1')
        with rinput.obresult.frames[0].open() as hdul:
            other = flow(hdul)
            f1 = other[0].data.copy() * iinfo[0]['texp'] * 1e-3

        self.logger.info('processing flat #2')
        with rinput.obresult.frames[1].open() as hdul:
            other = flow(hdul)
            f2 = other[0].data.copy() * iinfo[1]['texp'] * 1e-3

        # Preprocess...

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
            self.logger.debug('iter %d', niter)
            ratio, m, sigma = cosmetics(
                f1, f2, m, lowercut=lowercut, uppercut=uppercut)

            if self.intermediate_results:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    fits.writeto('numina-cosmetics-i%02d.fits' %
                                 niter, ratio, overwrite=True)
                    fits.writeto('numina-mask-i%02d.fits' % niter, m, overwrite=True)
                    fits.writeto('numina-sigma-i%02d.fits' %
                                 niter, m * 0.0 + sigma, overwrite=True)
            self.logger.debug(
                'iter %d, invalid points in input mask: %d', niter, ninvalid)
            self.logger.debug('iter %d, estimated sigma is %f', niter, sigma)
            n_ninvalid = numpy.count_nonzero(m)

            # Probably there is something wrong here
            # too much defective pixels
            if ninvalid / m.size >= 0.10:
                # This should set a flag in the output
                msg = 'defective pixels are greater than 10%'
                self.logger.warning(msg)

            if n_ninvalid == ninvalid:
                self.logger.info('convergence reached after %d iterations', niter)
                break
            self.logger.info('new invalid points: %d', n_ninvalid - ninvalid)
            ninvalid = n_ninvalid
        else:
            # This should set a flag in the output
            msg = 'convergence not reached after %d iterations' % maxiter
            self.logger.warning(msg)

        self.logger.info(
            'number of dead pixels %d', numpy.count_nonzero(m == PIXEL_DEAD))
        self.logger.info(
            'number of hot pixels %d', numpy.count_nonzero(m == PIXEL_HOT))

        if self.intermediate_results:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fits.writeto('numina-cosmetics.fits', ratio, overwrite=True)
                fits.writeto('numina-mask.fits', m, overwrite=True)
                fits.writeto(
                    'numina-sigma.fits', sigma * numpy.ones_like(m), overwrite=True)

        hdu = fits.PrimaryHDU(ratio)
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        ratiohdl = fits.HDUList([hdu])

        maskhdu = fits.PrimaryHDU(m)
        hdr = maskhdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        maskhdl = fits.HDUList([maskhdu])

        res = self.create_result(ratioframe=ratiohdl, maskframe=maskhdl)
        return res
