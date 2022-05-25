#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Beam switched-nodded image mode recipe of EMIR"""


from numina.core import Parameter
from numina.core import RecipeInput
from numina.core.recipeinout import define_input, define_result
from numina.core import Result, DataFrameType, RecipeResult
from numina.core.requirements import ObservationResultRequirement

from emirdrp.products import SourcesCatalog
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import Offsets_Requirement
from emirdrp.requirements import SkyImageSepTime_Requirement

from .shared import DirectImageCommon


class NBImageRecipeInput(RecipeInput):
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    extinction = Parameter(0.0, 'Mean atmospheric extinction')
    sources = Parameter(
        [], 'List of x, y coordinates to measure FWHM', optional=True)
    offsets = Offsets_Requirement()
    sky_images = Parameter(5, 'Images used to estimate the background'
                           ' before and after current image')
    sky_images_sep_time = SkyImageSepTime_Requirement()
    check_photometry_levels = Parameter(
        [0.5, 0.8], 'Levels to check the flux of the objects')
    check_photometry_actions = Parameter(
        ['warn', 'warn', 'default'], 'Actions to take on images')


class NBImageRecipeResult(RecipeResult):
    frame = Result(DataFrameType)
    catalog = Result(SourcesCatalog)


@define_input(NBImageRecipeInput)
@define_result(NBImageRecipeResult)
class NBImageRecipe(DirectImageCommon):

    """
    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing the TS in cycles
    between two, or more, sky positions. Displacements are larger
    than the EMIR FOV, so the images have no common area. Used
    for sky subtraction


    **Observing modes:**

        * Nodded/Beamswitched images

    """

    def run(self, ri):

        frame, catalog = self.process(ri,
                                      window=None, subpix=1,
                                      target_is_sky=False,
                                      stop_after=DirectImageCommon.FULLRED)

        return self.create_result(frame=frame, catalog=catalog)
