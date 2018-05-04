#
# Copyright 2016-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""
Spectroscopy mode, Stare Spectra
"""


from numina.core import Result
from numina.array.combine import median

from emirdrp.core.recipe import EmirRecipe
from emirdrp.processing.combine import basic_processing_with_combination
import emirdrp.requirements as reqs
import emirdrp.products as prods


class SkySpecRecipe(EmirRecipe):
    """Recipe to process data taken in spectral sky mode.

    """

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()

    skyspec = Result(prods.SkySpectrum)


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
