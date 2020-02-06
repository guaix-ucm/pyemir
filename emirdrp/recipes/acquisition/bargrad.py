#
# Copyright 2015-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Bar characterization using gradients for EMIR"""

from __future__ import division

from numina.core import Requirement, Result, Parameter
from numina.exceptions import RecipeError
import numina.types.array as tarray

import emirdrp.datamodel as datamodel
from emirdrp.processing.bars import slits_to_ds9_reg, find_bars
from emirdrp.core import EMIR_NBARS
from emirdrp.core.recipe import EmirRecipe
from numina.processing.combine import basic_processing_with_combination

import emirdrp.requirements as reqs
import emirdrp.products as prods


class BarDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement()

    bars_nominal_positions = Requirement(prods.NominalPositions,
                                         'Nominal positions of the bars'
                                         )
    median_filter_size = Parameter(5, 'Size of the median box')
    average_box_row_size = Parameter(7, 'Number of rows to average for fine centering (odd)')
    average_box_col_size = Parameter(21, 'Number of columns to extract for fine centering (odd)')
    fit_peak_npoints = Parameter(3, 'Number of points to use for fitting the peak (odd)')

    # Recipe Products
    frame = Result(prods.ProcessedImage)
    # derivative = Result(prods.ProcessedImage)
    slits = Result(tarray.ArrayType)
    positions3 = Result(tarray.ArrayType)
    positions5 = Result(tarray.ArrayType)
    positions7 = Result(tarray.ArrayType)
    positions9 = Result(tarray.ArrayType)
    DTU = Result(tarray.ArrayType)
    ROTANG = Result(float)
    TSUTC1 = Result(float)
    csupos = Result(tarray.ArrayType)
    csusens = Result(tarray.ArrayType)

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
            mecs_hdr = datamodel.get_mecs_header(hdulist)
            dtub, dtur = datamodel.get_dtur_from_header(mecs_hdr)
            csupos = datamodel.get_csup_from_header(mecs_hdr)
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
