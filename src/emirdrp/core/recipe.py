#
# Copyright 2008-2025 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import logging

import numina.core.recipes as recipes

import emirdrp.core.correctors as cor
import emirdrp.products as prods
import emirdrp.datamodel


class EmirRecipe(recipes.BaseRecipe):
    """Base class for all EMIR Recipes"""

    logger = logging.getLogger(__name__)
    datamodel = emirdrp.datamodel.EmirDataModel()

    def types_getter(self):
        imgtypes = [
            prods.MasterBadPixelMask,
            prods.MasterBias,
            prods.MasterDark,
            prods.MasterIntensityFlat,
            prods.MasterSpectralFlat,
            prods.MasterSky,
        ]
        getters = [
            cor.get_corrector_p,
            cor.get_corrector_b,
            cor.get_corrector_d,
            [cor.get_corrector_f, cor.get_checker],
            cor.get_corrector_sf,
            cor.get_corrector_s,
        ]

        return imgtypes, getters

    def create_getters(self):
        imgtypes, getters = self.types_getter()
        getters = self.get_filters(imgtypes, getters)
        getters_seq = [getters]
        return getters_seq

    def aggregate_result(self, result, rinput):
        return result

    def init_filters(self, rinput, ins=None):
        # TODO: we overload this method to
        # avoid rewriting all the calls like
        # flow = self.init_filters(rinput)
        # to
        # flow = self.init_filters(rinput)[0]
        # we should do it in the future
        return super().init_filters(rinput, ins)[0]
