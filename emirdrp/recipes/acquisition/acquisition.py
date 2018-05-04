#
# Copyright 2008-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""Recipe for the processing of target acquisition images."""


from emirdrp.core.recipe import EmirRecipe
import numina.core.dataholders as dh
import emirdrp.requirements as reqs
import emirdrp.products as prods


class TargetAcquisitionRecipe(EmirRecipe):

    """
    Acquire a target.

    Recipe for the processing of target acquisition images.

    **Observing modes:**

        * Target acquisition


    """

    # Requirements
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    
    # Products
    telescope_offset = dh.Result(prods.TelescopeOffset)

    def run(self, rinput):
        return self.create_result(telescope_offset=prods.TelescopeOffset())
