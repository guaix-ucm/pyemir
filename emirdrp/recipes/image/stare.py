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

import datetime

from numina.array.combine import median
from numina.core import Parameter
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement

from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.products import SourcesCatalog
import emirdrp.requirements as reqs
from emirdrp.processing.combine import basic_processing_with_combination
import emirdrp.processing.info as info
from .shared import DirectImageCommon


def timeit(method):
    def timed_method(self, rinput):

        time_start = datetime.datetime.utcnow()
        result = method(self, rinput)
        time_end = datetime.datetime.utcnow()
        result.time_it(time_start, time_end)
        return result

    return timed_method


class StareImageBaseRecipe(EmirRecipe):
    """Process images in Stare Image Mode"""

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement(optional=True)

    frame = Product(DataFrameType)

    @timeit
    def run(self, rinput):
        print(info.gather_info(rinput))
        self.logger.info('starting stare image reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow, method=median)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        # Update SEC to 0
        hdr['SEC'] = 0
        self.logger.info('end stare image reduction')
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
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    extinction = reqs.Extinction_Requirement()
    sources = reqs.Catalog_Requirement()
    offsets = reqs.Offsets_Requirement()
    iterations = Parameter(4, 'Iterations of the recipe')

    frame = Product(DataFrameType)
    catalog = Product(SourcesCatalog)

    def run(self, recipe_input):

        frame, catalog = self.process(recipe_input,
                                      window=None, subpix=1,
                                      stop_after=DirectImageCommon.PRERED)

        return self.create_result(frame=frame, catalog=catalog)
