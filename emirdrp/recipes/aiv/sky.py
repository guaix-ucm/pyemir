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

"""AIV Recipes for EMIR"""

import logging

from numina.array.combine import median
from numina.core import Product, Requirement
from numina.core.requirements import ObservationResultRequirement

from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType, MasterIntensityFlat
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.processing.combine import basic_processing_with_combination


_logger = logging.getLogger('numina.recipes.emir')


class TestSkyCorrectRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = Requirement(MasterIntensityFlat, 'Master Sky calibration')

    frame = Product(DataFrameType)

    def run(self, rinput):
        _logger.info('starting simple sky reduction')

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
