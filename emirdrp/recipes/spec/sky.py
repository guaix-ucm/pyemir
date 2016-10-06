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
Spectroscopy mode, Stare Spectra
"""


from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
from numina.array.combine import median
import emirdrp.requirements as reqs
from emirdrp.core import EmirRecipe
import emirdrp.products as prods
from emirdrp.processing.combine import basic_processing_with_combination

class SkySpecRecipe(EmirRecipe):
    """Recipe to process data taken in spectral sky mode.

    """

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()

    skyspec = Product(prods.SkySpectrum)


    def run(self, rinput):
        self.logger.info('starting spectral sky reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=True)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        self.logger.info('end sky spectral reduction')

        result = self.create_result(skyspec=hdulist)

        return result
