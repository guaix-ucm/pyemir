#
# Copyright 2008-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""Recipes for finding the best focus."""

from numina.core import Parameter
from numina.core import Result

import emirdrp.requirements as reqs
import emirdrp.products as prods
from emirdrp.core.recipe import EmirRecipe


class TelescopeRoughFocusRecipe(EmirRecipe):
    """Recipe to compute the telescope focus.

    **Observing modes:**

     * Telescope rough focus
     * Emir focus control

    **Inputs:**

     * A list of images
     * A list of sky images
     * Bias, dark, flat
     * A model of the detector
     * List of focii

    **Outputs:**
     * Best focus
    """

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    objects = Parameter([], 'List of x-y pair of object coordinates'),
    focus_range = Parameter([], 'Focus range: begin, end and step')

    focus = Result(prods.TelescopeFocus)

    def run(self, rinput):
        return self.create_result(focus=prods.TelescopeFocus())


class TelescopeFineFocusRecipe(EmirRecipe):
    """
    Recipe to compute the telescope focus.

    **Observing modes:**

        * Telescope fine focus

    **Inputs:**

     * A list of images
     * A list of sky images
     * Bias, dark, flat
     * A model of the detector
     * List of focii

    **Outputs:**
     * Best focus

    """

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    objects = Parameter([], 'List of x-y pair of object coordinates'),

    focus = Result(prods.TelescopeFocus)

    def run(self, rinput):
        return self.create_result(focus=prods.TelescopeFocus())


class DTUFocusRecipe(EmirRecipe):
    """
    Recipe to compute the DTU focus.

    **Observing modes:**

        * EMIR focus control

    **Inputs:**

     * A list of images
     * A list of sky images
     * Bias, dark, flat
     * A model of the detector
     * List of focii

    **Outputs:**
     * Best focus

    """

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    objects = Parameter([], 'List of x-y pair of object coordinates'),
    msm_pattern = Parameter([], 'List of x-y pair of slit coordinates'),
    dtu_focus_range = Parameter([], 'Focus range of the DTU: begin, end and step')

    focus = Result(prods.DTUFocus)


    def run(self, rinput):
        return self.create_result(focus=prods.DTUFocus())
