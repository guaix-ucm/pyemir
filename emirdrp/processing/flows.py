#
# Copyright 2013-2016 Universidad Complutense de Madrid
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

"""Load basic calibrations and create reduction flows"""

from __future__ import division
#
import logging

import numpy
from astropy.io import fits

from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import SkyCorrector
from numina.flow import SerialFlow
from numina.flow.node import IdNode
from numina.array import combine
from numina.array import combine_shape
#
import emirdrp.ext.gtc
from emirdrp.core import EMIR_BIAS_MODES
from emirdrp.core import gather_info
from emirdrp.core import offsets_from_wcs
from emirdrp.processing.badpixels import BadPixelCorrectorEmir
from emirdrp.processing.flatfield import FlatFieldCorrector
from emirdrp.processing.checkers import Checker


_logger = logging.getLogger('numina.recipes.emir')


def init_filters_generic(rinput, getters):
    # with BPM, bias, dark, flat and sky
    if emirdrp.ext.gtc.RUN_IN_GTC:
        _logger.debug('running in GTC environment')
    else:
        _logger.debug('running outside of GTC environment')

    meta = gather_info(rinput)
    _logger.debug('obresult info is %s', meta['obresult'])
    correctors = [getter(rinput, meta) for getter in getters]

    flow = SerialFlow(correctors)

    return flow


def get_corrector_p(rinput, meta):

    bpm_info = meta.get('master_bpm')
    if bpm_info is not None:
        with rinput.master_bpm.open() as hdul:
            _logger.info('loading BPM')
            _logger.debug('BPM image: %s', bpm_info)
            mbpm = hdul[0].data
            bpm_corrector = BadPixelCorrectorEmir(mbpm)
    else:
        _logger.info('BPM not provided, ignored')
        bpm_corrector = IdNode()

    return bpm_corrector


def get_corrector_b(rinput, meta):
    iinfo = meta['obresult']
    _logger.debug('images info: %s', iinfo)
    if iinfo:
        mode = iinfo[0]['readmode']
        if mode.lower() in EMIR_BIAS_MODES:
            use_bias = True
            _logger.info('readmode is %s, bias required', mode)
        else:
            use_bias = False
            _logger.info('readmode is %s, bias not required', mode)
    else:
        raise ValueError('cannot gather images info')

    # Loading calibrations
    if use_bias:
        bias_info = meta['master_bias']
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            _logger.debug('bias info: %s', bias_info)
            mbias = hdul[0].data
            bias_corrector = BiasCorrector(mbias)
    else:
        _logger.info('ignoring bias')
        bias_corrector = IdNode()

    return bias_corrector


def get_corrector_s(rinput, meta):
    sky_info = meta['master_sky']


    with rinput.master_sky.open() as msky_hdul:
        _logger.info('loading sky')
        _logger.debug('sky info: %s', sky_info)
        msky = msky_hdul[0].data
        sky_corrector = SkyCorrector(msky)

    return sky_corrector


def get_corrector_f(rinput, meta):
    flat_info = meta['master_flat']

    with rinput.master_flat.open() as mflat_hdul:
        _logger.info('loading intensity flat')
        _logger.debug('flat info: %s', flat_info)
        mflat = mflat_hdul[0].data
        # Check NaN and Ceros
        mask1 = mflat < 0
        mask2 = ~numpy.isfinite(mflat)
        if numpy.any(mask1):
            _logger.warning('flat has %d values below 0', mask1.sum())
        if numpy.any(mask2):
            _logger.warning('flat has %d NaN', mask2.sum())
        flat_corrector = FlatFieldCorrector(mflat)

    return flat_corrector


def get_corrector_d(rinput, meta):
    key = 'master_dark'
    CorrectorClass = DarkCorrector

    info = meta[key]
    req = getattr(rinput, key)
    with req.open() as hdul:
        _logger.info('loading %s', key)
        _logger.debug('%s info: %s', key, info)
        datac = hdul[0].data
        corrector = CorrectorClass(datac)

    return corrector


def get_checker(rinput, meta):

    return Checker()

