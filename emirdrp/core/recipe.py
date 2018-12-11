#
# Copyright 2008-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import logging

import numina.ext.gtc
import numina.core.recipes as recipes
import numina.core.recipeinout as recipeio
import numina.core.dataholders as dh
import numina.types.qc as qctype
import numina.types.obsresult as obtype
import numina.types.dataframe as dataframe
import numina.util.flow as flowmod

import emirdrp.core.correctors as cor
import emirdrp.processing.info
import emirdrp.products as prods
import emirdrp.datamodel


class EmirRecipeResult(recipeio.RecipeResult):
    """Specific RecipeResult for EMIR"""

    def time_it(self, time1, time2):
        values = self.attrs()
        for k, spec in self.stored().items():
            value = values[k]
            # Store for Images..
            if isinstance(value, dataframe.DataFrame):
                hdul = value.open()
                self.add_computation_time(hdul, time1, time2)

    def add_computation_time(self, img, time1, time2):
        img[0].header['NUMUTC1'] = time1.isoformat()
        img[0].header['NUMUTC2'] = time2.isoformat()
        return img


class EmirRecipe(recipes.BaseRecipe):
    """Base clase for all EMIR Recipes


    Attributes
    ----------
    qc : QualityControl, result, QC.GOOD by default

    logger :
         recipe logger

    datamodel : EmirDataModel

    """
    RecipeResult = EmirRecipeResult

    qc = dh.Result(obtype.QualityControlProduct,
                   destination='qc', default=qctype.QC.GOOD)
    logger = logging.getLogger('numina.recipes.emir')
    datamodel = emirdrp.datamodel.EmirDataModel()

    def types_getter(self):
        imgtypes = [prods.MasterBadPixelMask,
                    prods.MasterBias,
                    prods.MasterDark,
                    prods.MasterIntensityFlat,
                    prods.MasterSky
                    ]
        getters = [cor.get_corrector_p, cor.get_corrector_b, cor.get_corrector_d,
                   [cor.get_corrector_f, cor.get_checker], cor.get_corrector_s]

        return imgtypes, getters

    def get_filters(self):
        import collections
        imgtypes, getters = self.types_getter()
        used_getters = []
        for rtype, getter in zip(imgtypes, getters):
            self.logger.debug('get_filters, %s  %s', rtype, getter)
            if rtype is None:
                # Unconditional
                if isinstance(getter, collections.Iterable):
                    used_getters.extend(getter)
                else:
                    used_getters.append(getter)
            else:
                # Search
                for key, val in self.RecipeInput.stored().items():
                    if isinstance(val.type, rtype):
                        if isinstance(getter, collections.Iterable):
                            used_getters.extend(getter)
                        else:
                            used_getters.append(getter)
                        break
                else:
                    pass
        return used_getters

    def init_filters_generic(self, rinput, getters, ins):
        # with BPM, bias, dark, flat and sky
        if numina.ext.gtc.check_gtc():
            self.logger.debug('running in GTC environment')
        else:
            self.logger.debug('running outside of GTC environment')

        meta = emirdrp.processing.info.gather_info(rinput)
        self.logger.debug('obresult info')
        for entry in meta['obresult']:
            self.logger.debug('frame info is %s', entry)
        correctors = [getter(rinput, meta, ins, self.datamodel) for getter in getters]
        flow = flowmod.SerialFlow(correctors)
        return flow

    def init_filters(self, rinput, ins='EMIR'):
        getters = self.get_filters()
        return self.init_filters_generic(rinput, getters, ins)

    def aggregate_result(self, result, rinput):
        return result
