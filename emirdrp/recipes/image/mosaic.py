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

"""Mosaic image mode recipe of EMIR"""


from numina.core import Parameter
from numina.core import DataFrame
from numina.core import Product, DataFrameType
from numina.core.requirements import ObservationResultRequirement

from emirdrp.core import EmirRecipe
from emirdrp.products import SourcesCatalog


class MosaicRecipe(EmirRecipe):

    """
    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing to a number of
    sky positions, with separations of the order of the EMIR FOV.
    This command is designed to fully cover a given area on the
    sky, but can also be used to point to a number of sky positions
    on which acquisition is only made at the beginning. Supersky
    frame(s) can be built from the image series.

    **Observing modes:**

        * Mosaic images

    """

    obresult = ObservationResultRequirement()
    # FIXME: this parameter is optional
    sources = Parameter([], 'List of x, y coordinates to measure FWHM')

    frame = Product(DataFrameType)
    catalog = Product(SourcesCatalog)

    def run(self, ri):
        return self.create_result(frame=DataFrame(None),
                                  catalog=SourcesCatalog())

#
