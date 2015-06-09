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
    Recipes for finding the best focus.

'''

import logging

from numina.core import Parameter
from numina.core import DataProductRequirement
from numina.core import Product

from emirdrp.core import EmirRecipe
from emirdrp.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emirdrp.dataproducts import TelescopeFocus
from emirdrp.dataproducts import DTUFocus
from emirdrp.dataproducts import MasterIntensityFlat

__all__ = ['TelescopeRoughFocusRecipe',
           'TelescopeFineFocusRecipe',
           'DTUFocusRecipe',
           ]

_logger = logging.getLogger('numina.recipes.emir')


class TelescopeRoughFocusRecipe(EmirRecipe):

    '''Recipe to compute the telescope focus.

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
    '''

    master_bpm = DataProductRequirement(
        MasterBadPixelMask, 'Master bad pixel mask')
    master_bias = DataProductRequirement(MasterBias, 'Master bias image')
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    master_flat = DataProductRequirement(
        MasterIntensityFlat, 'Master intensity flatfield')
    objects = Parameter([], 'List of x-y pair of object coordinates'),
    focus_range = Parameter([], 'Focus range: begin, end and step')

    focus = Product(TelescopeFocus)

    def run(self, obresult, reqs):
        return self.create_result(focus=TelescopeFocus())


class TelescopeFineFocusRecipe(EmirRecipe):

    '''
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

    '''

    master_bpm = DataProductRequirement(
        MasterBadPixelMask, 'Master bad pixel mask')
    master_bias = DataProductRequirement(MasterBias, 'Master bias image')
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    master_flat = DataProductRequirement(
        MasterIntensityFlat, 'Master intensity flatfield')
    objects = Parameter([], 'List of x-y pair of object coordinates'),

    focus = Product(TelescopeFocus)

    def run(self, obresult, reqs):
        return self.create_result(focus=TelescopeFocus())


class DTUFocusRecipe(EmirRecipe):

    '''
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

    '''

    master_bpm = DataProductRequirement(
        MasterBadPixelMask, 'Master bad pixel mask')
    master_bias = DataProductRequirement(MasterBias, 'Master bias image')
    master_dark = DataProductRequirement(MasterDark, 'Master dark image')
    master_flat = DataProductRequirement(
        MasterIntensityFlat, 'Master intensity flatfield')
    objects = Parameter([], 'List of x-y pair of object coordinates'),
    msm_pattern = Parameter([], 'List of x-y pair of slit coordinates'),
    dtu_focus_range = Parameter(
        'dtu_focus_range', [], 'Focus range of the DTU: begin, end and step')

    focus = Product(DTUFocus)


    def run(self, obresult, reqs):
        return self.create_result(focus=DTUFocus())
