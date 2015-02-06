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

'''
Recipe for the processing of target acquisition images.

**Observing modes:**

    * Target acquisition
'''

import logging


from numina.core import Product

from numina.core.requirements import ObservationResultRequirement

from emir.core import EmirRecipe
from emir.dataproducts import TelescopeOffset, MSMPositions
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterBadPixelMaskRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement


__all__ = ['TargetAcquisitionRecipe', 'MaskImagingRecipe', 'MaskCheckRecipe']

_logger = logging.getLogger('emir.recipes')


class TargetAcquisitionRecipe(EmirRecipe):

    '''
    Acquire a target.

    Recipe for the processing of target acquisition images.

    **Observing modes:**

        * Target acquisition


    '''

    # Requirements
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    
    # Products
    telescope_offset = Product(TelescopeOffset)

    def run(self, obresult, reqs):
        return self.create_result(telescope_offset=TelescopeOffset())


class MaskImagingRecipe(EmirRecipe):

    '''Acquire a target.

    Mask image Recipe.

    Recipe to process mask images.

    **Observing modes:**

      *  Mask imaging
    '''

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()

    msm_positions = Product(MSMPositions)

    def run(self, obresult, reqs):
        return self.create_result(msm_positions=MSMPositions())


class MaskCheckRecipe(EmirRecipe):

    '''
    Acquire a target.

    Recipe for the processing of multi-slit/long-slit check images.

    **Observing modes:**

        * MSM and LSM check

    '''

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()

    msm_positions = Product(MSMPositions)
    telescope_offset = Product(TelescopeOffset)


    def run(self, obresult, reqs):
        return self.create_result(msm_positions=MSMPositions(),
                                  telescope_offset=TelescopeOffset())
