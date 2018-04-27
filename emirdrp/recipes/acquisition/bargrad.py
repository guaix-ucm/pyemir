#
# Copyright 2015-2018 Universidad Complutense de Madrid
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

"""Bar characterization using gradients for EMIR"""

from __future__ import division

from numina.core import Requirement, Product, Parameter
from numina.exceptions import RecipeError
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement

import emirdrp.datamodel as datamodel
from emirdrp.processing.bars import slits_to_ds9_reg, find_bars
from emirdrp.core import EMIR_NBARS
from emirdrp.core.recipe import EmirRecipe
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.products import DataFrameType, NominalPositions
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement


class BarDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    bars_nominal_positions = Requirement(NominalPositions,
                                         'Nominal positions of the bars'
                                         )
    median_filter_size = Parameter(5, 'Size of the median box')
    average_box_row_size = Parameter(7, 'Number of rows to average for fine centering (odd)')
    average_box_col_size = Parameter(21, 'Number of columns to extract for fine centering (odd)')
    fit_peak_npoints = Parameter(3, 'Number of points to use for fitting the peak (odd)')

    # Recipe Products
    frame = Product(DataFrameType)
    # derivative = Product(DataFrameType)
    slits = Product(ArrayType)
    positions3 = Product(ArrayType)
    positions5 = Product(ArrayType)
    positions7 = Product(ArrayType)
    positions9 = Product(ArrayType)
    DTU = Product(ArrayType)
    ROTANG = Product(float)
    TSUTC1 = Product(float)
    csupos = Product(ArrayType)
    csusens = Product(ArrayType)

    def run(self, rinput):
        self.logger.info('starting processing for bars detection')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(hdulist, 'reduced_image.fits')

        try:
            rotang = hdr['ROTANG']
            tsutc1 = hdr['TSUTC1']
            dtub, dtur = datamodel.get_dtur_from_header(hdr)
            csupos = datamodel.get_csup_from_header(hdr)
            if len(csupos) != 2 * EMIR_NBARS:
                raise RecipeError('Number of CSUPOS != 2 * NBARS')
            csusens = datamodel.get_cs_from_header(hdr)

        except KeyError as error:
            self.logger.error(error)
            raise RecipeError(error)

        self.logger.debug('start finding bars')
        allpos, slits = find_bars(hdulist,
                                  rinput.bars_nominal_positions,
                                  csupos,
                                  dtur,
                                  average_box_row_size=rinput.average_box_row_size,
                                  average_box_col_size=rinput.average_box_col_size,
                                  fit_peak_npoints=rinput.fit_peak_npoints,
                                  median_filter_size=rinput.median_filter_size,
                                  logger=self.logger
                                  )

        self.logger.debug('end finding bars')

        if self.intermediate_results:
            with open('ds9.reg', 'w') as ds9reg:
                slits_to_ds9_reg(ds9reg, slits)

        result = self.create_result(frame=hdulist,
                                    slits=slits,
                                    positions9=allpos[9],
                                    positions7=allpos[7],
                                    positions5=allpos[5],
                                    positions3=allpos[3],
                                    DTU=dtub,
                                    ROTANG=rotang,
                                    TSUTC1=tsutc1,
                                    csupos=csupos,
                                    csusens=csusens,
                                    )
        return result
