#
# Copyright 2008-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import logging
import collections.abc

import numina.ext.gtc
import numina.core.recipes as recipes
import numina.util.flow as flowmod

import emirdrp.core.correctors as cor
import emirdrp.processing.info
import emirdrp.products as prods
import emirdrp.datamodel


class EmirRecipe(recipes.BaseRecipe):
    """Base class for all EMIR Recipes"""
    logger = logging.getLogger(__name__)
    datamodel = emirdrp.datamodel.EmirDataModel()

    def types_getter(self):
        imgtypes = [prods.MasterBadPixelMask,
                    prods.MasterBias,
                    prods.MasterDark,
                    prods.MasterIntensityFlat,
                    prods.MasterSpectralFlat,
                    prods.MasterSky
                    ]
        getters = [cor.get_corrector_p, cor.get_corrector_b, cor.get_corrector_d,
                   [cor.get_corrector_f, cor.get_checker], cor.get_corrector_sf, cor.get_corrector_s]

        return imgtypes, getters

    def get_filters(self, imgtypes, getters):

        used_getters = []
        for rtype, getter in zip(imgtypes, getters):
            self.logger.debug('get_filters, %s  %s', rtype, getter)
            if rtype is None:
                # Unconditional
                if isinstance(getter, collections.abc.Iterable):
                    used_getters.extend(getter)
                else:
                    used_getters.append(getter)
            else:
                # Search
                for key, val in self.RecipeInput.stored().items():
                    if isinstance(val.type, rtype):
                        if isinstance(getter, collections.abc.Iterable):
                            used_getters.extend(getter)
                        else:
                            used_getters.append(getter)
                        break
                else:
                    pass
        return used_getters

    def init_filters_generic(self, rinput, getters_seq, ins):
        # with BPM, bias, dark, flat and sky
        if numina.ext.gtc.check_gtc():
            self.logger.debug('running in GTC environment')
        else:
            self.logger.debug('running outside of GTC environment')

        meta = self.gather_info(rinput)
        self.logger.debug('obresult info')
        for entry in meta['obresult']:
            self.logger.debug('frame info is %s', entry)

        reduction_flows = []
        for getters in getters_seq:
            correctors = [getter(rinput, meta, ins, self.datamodel)
                          for getter in getters]
            flow = flowmod.SerialFlow(correctors)
            reduction_flows.append(flow)

        return reduction_flows

    def init_filters(self, rinput, ins=None):
        if ins is None:
            ins = rinput.obresult.configuration
        imgtypes, getters = self.types_getter()
        getters = self.get_filters(imgtypes, getters)
        getters_seq = [getters]
        # FIXME:
        return self.init_filters_generic(rinput, getters_seq, ins)[0]

    def aggregate_result(self, result, rinput):
        return result

    def gather_info(self, rinput):
        meta = emirdrp.processing.info.gather_info(rinput)
        return meta
