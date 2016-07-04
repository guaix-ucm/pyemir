#
# Copyright 2011-2014 Universidad Complutense de Madrid
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
Image mode recipes of EMIR
"""

import logging

from numina.array.combine import median
from numina.core import Parameter
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement

from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.products import SourcesCatalog
from emirdrp.requirements import Catalog_Requirement
from emirdrp.requirements import Extinction_Requirement
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement
from emirdrp.requirements import Offsets_Requirement
from emirdrp.processing.combine import basic_processing_with_combination
from .shared import DirectImageCommon

_logger = logging.getLogger('numina.recipes.emir')


class StareImageBaseRecipe(EmirRecipe):
    """Process images in Stare Image Mode"""

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    frame = Product(DataFrameType)

    def run(self, rinput):
        _logger.info('starting stare image reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow, method=median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        # Update SEC to 0
        hdr['SEC'] = 0
        _logger.info('end stare image reduction')
        result = self.create_result(frame=hdulist)

        return result


class StareImageRecipe(DirectImageCommon):

    """
    The effect of recording images of the sky in a given pointing
    position of the TS


    **Observing modes:**

        * Stare image

    """

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    extinction = Extinction_Requirement()
    sources = Catalog_Requirement()
    offsets = Offsets_Requirement()
    iterations = Parameter(4, 'Iterations of the recipe')

    frame = Product(DataFrameType)
    catalog = Product(SourcesCatalog)

    def run(self, recipe_input):

        frame, catalog = self.process(recipe_input,
                                      window=None, subpix=1,
                                      stop_after=DirectImageCommon.PRERED)

        return self.create_result(frame=frame, catalog=catalog)
