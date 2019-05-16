#
# Copyright 2013-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""AIV Recipes for EMIR"""


from numina.array.combine import median
from numina.core import Result, Requirement
from numina.processing.combine import basic_processing_with_combination

from emirdrp.core.recipe import EmirRecipe
import emirdrp.requirements as reqs
import emirdrp.products as prods


class TestSkyCorrectRecipe(EmirRecipe):

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = Requirement(prods.MasterIntensityFlat, 'Master Sky calibration')

    frame = Result(prods.ProcessedImage)

    def run(self, rinput):
        self.logger.info('starting simple sky reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow, method=median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        # Update SEC to 0
        hdr['SEC'] = 0
        
        result = self.create_result(frame=hdulist)

        return result

from numina.core import ObservationResult


class StareImageRecipeInputBuilder(object):
    '''Class to build StareImageRecipe inputs from the Observation Results.

       Fetches SKY calibration image from the archive


    '''

    def __init__(self, dal):
        self.dal = dal
        self.sky_image = None

    def buildRecipeInput(self, obsres):

        if self.sky_image is None:
            print('obtaining SKY image')
            sky_cal_result = self.dal.getLastRecipeResult("EMIR", "EMIR", "IMAGE_SKY")
            self.sky_image = sky_cal_result['elements']['skyframe']

        obsres['master_sky'] = self.sky_image
        newOR = ObservationResult()
        newOR.frames = obsres['frames']
        obsres['obresult'] = newOR
        newRI = StareImageRecipeInput(**obsres)

        return newRI


StareImageRecipeInput = TestSkyCorrectRecipe.RecipeInput
TestSkyCorrectRecipe.InputBuilder = StareImageRecipeInputBuilder
