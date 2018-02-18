#
# Copyright 2016 Universidad Complutense de Madrid
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

"""
Spectroscopy mode, coadd ABBA
"""


import numina.core
import numina.array.combine as combine
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
import numina.exceptions

from emirdrp.core import EmirRecipe
import emirdrp.products as prods
from emirdrp.processing.combine import basic_processing_with_combination


class CoaddABBARecipe(EmirRecipe):
    """Process images in ABBA mode"""

    obresult = ObservationResultRequirement()

    spec_coadd_abba = Product(prods.DataFrameType)

    @classmethod
    def build_recipe_input(cls, obsres, dal, pipeline='default'):
        return cls.build_recipe_input_gtc(obsres, dal, pipeline=pipeline)

    @classmethod
    def build_recipe_input_gtc(cls, obsres, dal, pipeline='default'):
        cls.logger.debug('start recipe input builder')
        print(dir(obsres))
        stareImagesIds = obsres.stareSpectraIds
        cls.logger.debug('ABBA images IDS %s: ', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['spec_abba'])

        newOR = numina.core.ObservationResult()
        newOR.frames = stareImages
        newRI = cls.create_input(obresult=newOR)
        cls.logger.debug('end recipe input builder')
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

        result = self.create_result(spec_coadd_abba=hdulist)
        self.logger.info('end spectroscopy ABBA coadd reduction')
        return result


class CoaddRecipe(EmirRecipe):
    """Generic Coadd Recipe"""

    obresult = ObservationResultRequirement()

    result_coadd = Product(prods.DataFrameType)

    @classmethod
    def build_recipe_input(cls, obsres, dal, pipeline='default'):
        return cls.build_recipe_input_gtc(obsres, dal, pipeline=pipeline)

    @classmethod
    def build_recipe_input_gtc(cls, obsres, dal, pipeline='default'):
        cls.logger.debug('start recipe input builder')

        # This depends on the RecipeResult
        result_field = 'spec_abba'

        stareImagesIds = obsres.stareSpectraIds
        cls.logger.debug('Coadd images IDS %s: ', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements'][result_field])

        newOR = numina.core.ObservationResult()
        newOR.frames = stareImages
        newRI = cls.create_input(obresult=newOR)
        cls.logger.debug('end recipe input builder')
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
