#
# Copyright 2013-2014 Universidad Complutense de Madrid
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

'''Recipe to detect slits in the AIV mask'''

from __future__ import division
#
import logging
import math

import numpy
import scipy.interpolate as itpl
import scipy.optimize as opz
from astropy.modeling import models, fitting
from astropy.io import fits

from numina import __version__
from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import FlatFieldCorrector, SkyCorrector
from numina.flow import SerialFlow
from numina.flow.node import IdNode
from numina.array import combine
#
from emir.core import EMIR_BIAS_MODES
from emir.core import gather_info


_logger = logging.getLogger('numina.recipes.emir')


def init_filters_bdfs(rinput):
    # with bias, dark, flat and sky
    meta = gather_info(rinput)
    iinfo = meta['obresult']

    if iinfo:
        mode = iinfo[0]['readmode']
        if mode.lower() in EMIR_BIAS_MODES:
            use_bias = True
            _logger.info('readmode is %s, bias required', mode)
        else:
            use_bias = False
            _logger.info('readmode is %s, bias not required', mode)
    else:
        raise ValueError('cannot gather frames info')

    dark_info = meta['master_dark']
    flat_info = meta['master_flat']
    sky_info = meta['master_sky']

    print('images info:', iinfo)
    if use_bias:
        bias_info = meta['master_bias']
        print('bias info:', bias_info)
        _logger.debug('bias info: %s', bias_info)

    print('dark info:', dark_info)
    _logger.debug('dark info: %s', dark_info)
    print('flat info:', flat_info)
    _logger.debug('flat info: %s', flat_info)
    print('sky info:', sky_info)
    _logger.debug('sky info: %s', sky_info)

    # Loading calibrations
    if use_bias:
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            mbias = hdul[0].data
            bias_corrector = BiasCorrector(mbias)
    else:
        _logger.info('ignoring bias')
        bias_corrector = IdNode()

    with rinput.master_dark.open() as mdark_hdul:
        _logger.info('loading dark')
        mdark = mdark_hdul[0].data
        dark_corrector = DarkCorrector(mdark)

    with rinput.master_flat.open() as mflat_hdul:
        _logger.info('loading intensity flat')
        mflat = mflat_hdul[0].data
        flat_corrector = FlatFieldCorrector(mflat)

    with rinput.master_sky.open() as msky_hdul:
        _logger.info('loading sky')
        msky = msky_hdul[0].data
        sky_corrector = SkyCorrector(msky)

    flow = SerialFlow([bias_corrector,
                       dark_corrector,
                       flat_corrector,
                       sky_corrector]
                      )

    return flow


def init_filters_bdf(rinput):
    # with bias, dark, flat and sky
    meta = gather_info(rinput)
    iinfo = meta['obresult']

    if iinfo:
        mode = iinfo[0]['readmode']
        if mode.lower() in EMIR_BIAS_MODES:
            use_bias = True
            _logger.info('readmode is %s, bias required', mode)
        else:
            use_bias = False
            _logger.info('readmode is %s, bias not required', mode)
    else:
        raise ValueError('cannot gather frames info')

    dark_info = meta['master_dark']
    flat_info = meta['master_flat']

    print('images info:', iinfo)
    if use_bias:
        bias_info = meta['master_bias']
        print('bias info:', bias_info)
        _logger.debug('bias info: %s', bias_info)

    print('dark info:', dark_info)
    _logger.debug('dark info: %s', dark_info)
    print('flat info:', flat_info)
    _logger.debug('flat info: %s', flat_info)

    # Loading calibrations
    if use_bias:
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            mbias = hdul[0].data
            bias_corrector = BiasCorrector(mbias)
    else:
        _logger.info('ignoring bias')
        bias_corrector = IdNode()

    with rinput.master_dark.open() as mdark_hdul:
        _logger.info('loading dark')
        mdark = mdark_hdul[0].data
        dark_corrector = DarkCorrector(mdark)

    with rinput.master_flat.open() as mflat_hdul:
        _logger.info('loading intensity flat')
        mflat = mflat_hdul[0].data
        flat_corrector = FlatFieldCorrector(mflat)

    flow = SerialFlow([bias_corrector,
                       dark_corrector,
                       flat_corrector
                       ]
                      )

    return flow


def init_filters_bd(rinput):
    # with bias, dark, flat and sky
    meta = gather_info(rinput)
    iinfo = meta['obresult']

    if iinfo:
        mode = iinfo[0]['readmode']
        if mode.lower() in EMIR_BIAS_MODES:
            use_bias = True
            _logger.info('readmode is %s, bias required', mode)
        else:
            use_bias = False
            _logger.info('readmode is %s, bias not required', mode)
    else:
        raise ValueError('cannot gather frames info')

    dark_info = meta['master_dark']

    print('images info:', iinfo)
    if use_bias:
        bias_info = meta['master_bias']
        print('bias info:', bias_info)
        _logger.debug('bias info: %s', bias_info)

    print('dark info:', dark_info)
    _logger.debug('dark info: %s', dark_info)

    # Loading calibrations
    if use_bias:
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            mbias = hdul[0].data
            bias_corrector = BiasCorrector(mbias)
    else:
        _logger.info('ignoring bias')
        bias_corrector = IdNode()

    with rinput.master_dark.open() as mdark_hdul:
        _logger.info('loading dark')
        mdark = mdark_hdul[0].data
        dark_corrector = DarkCorrector(mdark)

    flow = SerialFlow([bias_corrector,
                       dark_corrector
                       ]
                      )

    return flow


def init_filters_b(rinput):
    # with bias, dark, flat and sky
    meta = gather_info(rinput)
    iinfo = meta['obresult']

    if iinfo:
        mode = iinfo[0]['readmode']
        if mode.lower() in EMIR_BIAS_MODES:
            use_bias = True
            _logger.info('readmode is %s, bias required', mode)
        else:
            use_bias = False
            _logger.info('readmode is %s, bias not required', mode)
    else:
        raise ValueError('cannot gather frames info')

    print('images info:', iinfo)
    if use_bias:
        bias_info = meta['master_bias']
        print('bias info:', bias_info)
        _logger.debug('bias info: %s', bias_info)

    # Loading calibrations
    if use_bias:
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            mbias = hdul[0].data
            bias_corrector = BiasCorrector(mbias)
    else:
        _logger.info('ignoring bias')
        bias_corrector = IdNode()

    flow = SerialFlow([bias_corrector])

    return flow


def basic_processing_with_combination(rinput, flow,
                                      method=combine.mean,
                                      errors=True):

    odata = []
    cdata = []
    try:
        _logger.info('processing input frames')
        for frame in rinput.obresult.frames:
            hdulist = frame.open()
            fname = hdulist.filename()
            if fname:
                _logger.info('input is %s', fname)
            else:
                _logger.info('input is %s', hdulist)

            final = flow(hdulist)
            _logger.debug('output is input: %s', final is hdulist)

            cdata.append(final)

            # Files to be closed at the end
            odata.append(hdulist)
            if final is not hdulist:
                odata.append(final)

        _logger.info("stacking %d images using 'mean'", len(cdata))
        data = method([d[0].data for d in cdata], dtype='float32')
        hdu = fits.PrimaryHDU(data[0], header=cdata[0][0].header.copy())
        hdu.header['NUMXVER'] = (__version__, 'Numina package version')
        if errors:
            varhdu = fits.ImageHDU(data[1], name='VARIANCE')
            num = fits.ImageHDU(data[2], name='MAP')
            result = fits.HDUList([hdu, varhdu, num])
        else:
            result = fits.HDUList([hdu])
    finally:
        _logger.debug('closing images')
        for hdulist in odata:
            hdulist.close()

    _logger.debug('update result header')

    return result
