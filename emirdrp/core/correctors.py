#
# Copyright 2008-2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


import logging

import numpy
import numina.util.node
import numina.datamodel as dm
import numina.processing as proc

import emirdrp.core as c

_logger = logging.getLogger(__name__)


def get_corrector_p(rinput, meta, ins, datamodel):
    key = 'master_bpm'
    info = meta.get(key)
    corrector_class = proc.BadPixelCorrector

    if info is not None:
        inputval = getattr(rinput, key)
        with inputval.open() as hdul:
            _logger.info('loading "%s"', key)
            _logger.debug('info: %s', info)
            corrector = corrector_class(
                hdul[0].data,
                datamodel=datamodel,
                calibid=dm.get_imgid(hdul, prefix=True)
            )
    else:
        _logger.info('"%s" not provided, ignored', key)
        corrector = numina.util.node.IdNode()

    return corrector


def get_corrector_b(rinput, meta, ins, datamodel):
    iinfo = meta['obresult']
    if iinfo:
        mode = iinfo[0]['readmode']
        if mode.lower() in c.EMIR_BIAS_MODES:
            use_bias = True
            _logger.info('readmode is %s, bias required', mode)
        else:
            use_bias = False
            _logger.info('readmode is %s, bias not required', mode)
    else:
        # raise ValueError('cannot gather images info')
        use_bias = False

    # Loading calibrations
    if use_bias:
        bias_info = meta['master_bias']
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            _logger.debug('bias info: %s', bias_info)
            mbias = hdul[0].data
            bias_corrector = proc.BiasCorrector(
                mbias,
                datamodel=datamodel,
                calibid=dm.get_imgid(hdul, prefix=True)
            )
    else:
        _logger.info('ignoring bias')
        bias_corrector = numina.util.node.IdNode()

    return bias_corrector


def get_corrector_s(rinput, meta, ins, datamodel):
    sky_info = meta.get('master_sky')

    if sky_info is None:
        return numina.util.node.IdNode()
    else:
        with rinput.master_sky.open() as hdul:
            _logger.info('loading sky')
            _logger.debug('sky info: %s', sky_info)
            sky_corrector = proc.SkyCorrector(
                hdul[0].data,
                datamodel=datamodel,
                calibid=dm.get_imgid(hdul, prefix=True)
            )
        return sky_corrector


def get_corrector_f(rinput, meta, ins, datamodel):
    """Corrector for intensity flat"""
    from emirdrp.processing.flatfield import FlatFieldCorrector
    flat_info = meta['master_flat']
    with rinput.master_flat.open() as hdul:
        _logger.info('loading intensity flat')
        _logger.debug('flat info: %s', flat_info)
        mflat = hdul[0].data
        # Check NaN and Ceros
        mask1 = mflat < 0
        mask2 = ~numpy.isfinite(mflat)
        if numpy.any(mask1):
            _logger.warning('flat has %d values below 0', mask1.sum())
        if numpy.any(mask2):
            _logger.warning('flat has %d NaN', mask2.sum())
        flat_corrector = FlatFieldCorrector(
            mflat,
            datamodel=datamodel,
            calibid=dm.get_imgid(hdul, prefix=True)
        )

    return flat_corrector


def get_corrector_sf(rinput, meta, ins, datamodel):
    """Corrector for spectral flat"""
    from emirdrp.processing.flatfield import FlatFieldCorrector
    flat_info = meta['master_flat']
    with rinput.master_flat.open() as hdul:
        _logger.info('loading spectral flat')
        _logger.debug('flat info: %s', flat_info)
        mflat = hdul[0].data
        # Check NaN and Ceros
        mask1 = mflat < 0
        mask2 = ~numpy.isfinite(mflat)
        if numpy.any(mask1):
            _logger.warning('flat has %d values below 0', mask1.sum())
        if numpy.any(mask2):
            _logger.warning('flat has %d NaN', mask2.sum())
        flat_corrector = FlatFieldCorrector(mflat,
                                            datamodel=datamodel,
                                            calibid=dm.get_imgid(hdul, prefix=True))

    return flat_corrector


def get_corrector_d(rinput, meta, ins, datamodel):
    key = 'master_dark'

    corrector = get_corrector_gen(rinput, datamodel, proc.DarkCorrector, key)
    return corrector


def get_corrector_gen(rinput, datamodel, CorrectorClass, key):
    value = getattr(rinput, key)
    req = rinput.stored()[key]
    if value is None:
        if req.optional:
            _logger.info('"%s" not provided, ignored', key)
            corrector = numina.util.node.IdNode()
            return corrector
        else:
            msg = '"{}" not provided, is required'.format(key)
            raise ValueError(msg)
    else:
        with value.open() as hdul:
            datac = hdul['primary'].data
            corrector = CorrectorClass(
                datac,
                calibid=dm.get_imgid(hdul, prefix=True),
                datamodel=datamodel
            )
        return corrector


def get_checker(rinput, meta, ins, datamodel):
    from emirdrp.processing.checkers import Checker
    return Checker()
