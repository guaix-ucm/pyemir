#
# Copyright 2008-2014 Universidad Complutense de Madrid
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

import collections
import logging
import numpy

from numina.flow.node import IdNode
from numina.flow.processing import BadPixelCorrector
from numina.core import BaseRecipe, Product, RecipeResult, DataFrame
from numina.core.products import QualityControlProduct

import emirdrp.products as prods
import emirdrp.ext.gtc
import emirdrp.processing.info
from emirdrp.processing.datamodel import EmirDataModel


_logger = logging.getLogger('numina.recipes.emir')


def get_corrector_p(rinput, meta):
    bpm_info = meta.get('master_bpm')
    if bpm_info is not None:
        with rinput.master_bpm.open() as hdul:
            _logger.info('loading BPM')
            _logger.debug('BPM image: %s', bpm_info)
            mbpm = hdul[0].data
            bpm_corrector = BadPixelCorrector(mbpm, datamodel=EmirDataModel())
    else:
        _logger.info('BPM not provided, ignored')
        bpm_corrector = IdNode()

    return bpm_corrector


def get_corrector_b(rinput, meta):
    from numina.flow.processing import BiasCorrector
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
        raise ValueError('cannot gather images info')

    # Loading calibrations
    if use_bias:
        bias_info = meta['master_bias']
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            _logger.debug('bias info: %s', bias_info)
            mbias = hdul[0].data
            bias_corrector = BiasCorrector(mbias, datamodel=EmirDataModel())
    else:
        _logger.info('ignoring bias')
        bias_corrector = IdNode()

    return bias_corrector


def get_corrector_s(rinput, meta):
    from numina.flow.processing import SkyCorrector
    sky_info = meta.get('master_sky')

    if sky_info is None:
        return IdNode()
    else:
        with rinput.master_sky.open() as msky_hdul:
            _logger.info('loading sky')
            _logger.debug('sky info: %s', sky_info)
            msky = msky_hdul[0].data
            sky_corrector = SkyCorrector(msky, datamodel=EmirDataModel())

        return sky_corrector


def get_corrector_f(rinput, meta):
    from emirdrp.processing.flatfield import FlatFieldCorrector
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
        flat_corrector = FlatFieldCorrector(mflat, datamodel=EmirDataModel())

    return flat_corrector


def get_corrector_d(rinput, meta):
    from numina.flow.processing import DarkCorrector
    key = 'master_dark'
    CorrectorClass = DarkCorrector
    datamodel = EmirDataModel()

    info = meta[key]
    _logger.debug('datamodel is %s', datamodel)
    req = getattr(rinput, key)
    with req.open() as hdul:
        _logger.info('loading %s', key)
        _logger.debug('%s info: %s', key, info)
        datac = hdul['primary'].data
        corrector = CorrectorClass(datac,
                                   calibid=datamodel.get_imgid(hdul),
                                   datamodel=datamodel)

    return corrector


def get_checker(rinput, meta):
    from emirdrp.processing.checkers import Checker
    return Checker()


class EmirRecipeResult(RecipeResult):

    def time_it(self, time1, time2):
        values = self.attrs()
        for k, spec in self.stored().items():
            value = values[k]
            # Store for Images..
            if isinstance(value, DataFrame):
                import emirdrp.processing.datamodel
                d = emirdrp.processing.datamodel.EmirDataModel()
                hdul = value.open()
                d.add_computation_time(hdul, time1, time2)


class EmirRecipe(BaseRecipe):

    RecipeResult = EmirRecipeResult

    qc = Product(QualityControlProduct, dest='qc')

    logger = logging.getLogger('numina.recipes.emir')

    def __init__(self, version="1"):
        super(EmirRecipe, self).__init__(version=version)

    @classmethod
    def types_getter(cls):
        imgtypes = [prods.MasterBadPixelMask,
                    prods.MasterBias,
                    prods.MasterDark,
                    prods.MasterIntensityFlat,
                    prods.MasterSky
                    ]
        getters = [get_corrector_p, get_corrector_b, get_corrector_d,
                   [get_corrector_f, get_checker], get_corrector_s]

        return imgtypes, getters

    @classmethod
    def load_getters(cls):
        imgtypes, getters = cls.types_getter()
        used_getters = []
        for rtype, getter in zip(imgtypes, getters):
            for key, val in cls.RecipeInput.stored().items():
                if isinstance(val.type, rtype):
                    if isinstance(getter, collections.Iterable):
                        used_getters.extend(getter)
                    else:
                        used_getters.append(getter)
                    break
            else:
                pass
        return used_getters

    @classmethod
    def init_filters_generic(cls, rinput, getters):
        from numina.flow import SerialFlow
        # with BPM, bias, dark, flat and sky
        if emirdrp.ext.gtc.RUN_IN_GTC:
            cls.logger.debug('running in GTC environment')
        else:
            cls.logger.debug('running outside of GTC environment')

        meta = emirdrp.processing.info.gather_info(rinput)
        cls.logger.debug('obresult info')
        for entry in meta['obresult']:
            cls.logger.debug('frame info is %s', entry)
        correctors = [getter(rinput, meta) for getter in getters]

        flow = SerialFlow(correctors)

        return flow

    @classmethod
    def init_filters(cls, rinput):
        getters = cls.load_getters()
        return cls.init_filters_generic(rinput, getters)


EMIR_BIAS_MODES = ['simple', 'bias', 'single']
EMIR_READ_MODES = ['simple', 'bias', 'single', 'cds', 'fowler', 'ramp']

EMIR_PIXSCALE = 18.0
EMIR_GAIN = 5.0 # ADU / e-
EMIR_RON = 5.69 # ADU
EMIR_NBARS = 55
