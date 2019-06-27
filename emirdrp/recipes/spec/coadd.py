#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Spectroscopy mode, coadd ABBA
"""


import numina.core
import numina.array.combine as combine
from numina.core import Result
from numina.core.requirements import ObservationResultRequirement
import numina.exceptions
import numina.ext.gtc
import numina.core.query as qmod
from numina.processing.combine import basic_processing_with_combination

from emirdrp.core.recipe import EmirRecipe
import emirdrp.products as prods


class CoaddABBARecipe(EmirRecipe):
    """Process images in ABBA mode"""

    obresult = ObservationResultRequirement(
        query_opts=qmod.ResultOf(
            'LS_ABBA.reduced_mos_abba',
            node='children',
            id_field="stareSpectraIds"
        )
    )

    reduced_mos_abba = Result(prods.ProcessedMOS)

    def build_recipe_input(self, obsres, dal):
        if numina.ext.gtc.check_gtc():
            self.logger.debug('running in GTC environment')
            return self.build_recipe_input_gtc(obsres, dal)
        else:
            self.logger.debug('running outside of GTC environment')
            return super(CoaddABBARecipe, self).build_recipe_input(
                obsres, dal
            )

    def build_recipe_input_gtc(self, obsres, dal):
        self.logger.debug('start recipe input builder')
        stareImagesIds = obsres.stareSpectraIds
        self.logger.debug('ABBA images IDS %s: ', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['reduced_mos_abba'])

        newOR = numina.core.ObservationResult()
        newOR.frames = stareImages
        newRI = self.create_input(obresult=newOR)
        self.logger.debug('end recipe input builder')
        return newRI

    def run(self, rinput):
        self.logger.info('starting spectroscopy ABBA coadd reduction')

        flow = self.init_filters(rinput)

        nimages = len(rinput.obresult.frames)
        self.logger.info('we receive %d images', nimages)
        if nimages == 0:
            msg = 'Received %d images' % nimages
            raise numina.exceptions.RecipeError(msg)

        hdulist = basic_processing_with_combination(
            rinput,
            flow,
            method=combine.mean,
            prolog="Process Coadd ABBA"
        )

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        # Update SEC to 0
        # hdr['SEC'] = 0

        result = self.create_result(reduced_mos_abba=hdulist)
        self.logger.info('end spectroscopy ABBA coadd reduction')
        return result


class CoaddRecipe(EmirRecipe):
    """Generic Coadd Recipe"""

    obresult = ObservationResultRequirement()

    result_coadd = Result(prods.ProcessedMOS)

    def build_recipe_input(self, obsres, dal):
        if numina.ext.gtc.check_gtc():
            self.logger.debug('running in GTC environment')
            return self.build_recipe_input_gtc(obsres, dal)
        else:
            self.logger.debug('running outside of GTC environment')
            return super(CoaddRecipe, self).build_recipe_input(
                obsres, dal
            )

    def build_recipe_input_gtc(self, obsres, dal):
        self.logger.debug('start recipe input builder')

        # This depends on the RecipeResult
        result_field = 'reduced_mos_abba'

        stareImagesIds = obsres.stareSpectraIds
        self.logger.debug('Coadd images IDS %s: ', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements'][result_field])

        newOR = numina.core.ObservationResult()
        newOR.frames = stareImages
        newRI = self.create_input(obresult=newOR)
        self.logger.debug('end recipe input builder')
        return newRI

    def run(self, rinput):
        self.logger.info('starting coadd reduction')

        flow = self.init_filters(rinput)

        nimages = len(rinput.obresult.frames)
        self.logger.info('we receive %d images', nimages)
        if nimages == 0:
            msg = 'Received %d images' % nimages
            raise numina.exceptions.RecipeError(msg)

        hdulist = basic_processing_with_combination(
            rinput,
            flow,
            method=combine.mean,
            prolog="Process Generic Coadd"
        )

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        # Update SEC to 0
        # hdr['SEC'] = 0

        result = self.create_result(spec_coadd_abba=hdulist)
        self.logger.info('end coadd reduction')
        return result
