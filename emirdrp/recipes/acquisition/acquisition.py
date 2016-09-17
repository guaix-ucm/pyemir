#
# Copyright 2008-2016 Universidad Complutense de Madrid
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
Recipe for the processing of target acquisition images.

**Observing modes:**

    * Target acquisition
"""

import logging


from numina.core import Product

from numina.core.requirements import ObservationResultRequirement

from emirdrp.core import EmirRecipe
from emirdrp.products import TelescopeOffset, MSMPositions
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement


_logger = logging.getLogger('emirdrp.recipes')


class TargetAcquisitionRecipe(EmirRecipe):

    """
    Acquire a target.

    Recipe for the processing of target acquisition images.

    **Observing modes:**

        * Target acquisition


    """

    # Requirements
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    
    # Products
    telescope_offset = Product(TelescopeOffset)

    def run(self, rinput):
        return self.create_result(telescope_offset=TelescopeOffset())
