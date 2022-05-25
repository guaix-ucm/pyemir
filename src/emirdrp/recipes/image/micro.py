#
# Copyright 2012-2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Micro-dithering Recipe of EMIR

"""

from numina.core import Parameter
from numina.core.recipeinout import define_input, define_result
from numina.core import Result, DataFrameType, RecipeInput, RecipeResult
from numina.core.requirements import ObservationResultRequirement

from emirdrp.products import SourcesCatalog
from emirdrp.requirements import Offsets_Requirement
from emirdrp.requirements import SkyImageSepTime_Requirement
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement

from .shared import DirectImageCommon


class MicroditheredImageRecipeInput(RecipeInput):
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    extinction = Parameter(0.0, 'Mean atmospheric extinction')
    sources = Parameter([],
                        'List of x, y coordinates to measure FWHM',
                        optional=True)
    offsets = Offsets_Requirement()
    iterations = Parameter(4, 'Iterations of the recipe')
    sky_images = Parameter(5, 'Images used to estimate the background before '
                           'and after current image')
    sky_images_sep_time = SkyImageSepTime_Requirement()
    check_photometry_levels = Parameter(
        [0.5, 0.8], 'Levels to check the flux of the objects')
    check_photometry_actions = Parameter(
        ['warn', 'warn', 'default'], 'Actions to take on images')
    subpixelization = Parameter(4, 'Number of subdivisions in each pixel side')
    window = Parameter([], 'Region of interesting data', optional=True)


class MicroditheredImageRecipeResult(RecipeResult):
    frame = Result(DataFrameType)
    catalog = Result(SourcesCatalog)


@define_input(MicroditheredImageRecipeInput)
@define_result(MicroditheredImageRecipeResult)
class MicroditheredImageRecipe(DirectImageCommon):

    """
    Recipe for the reduction of microdithering imaging.

    The effect of recording a series of stare images, with the same
    acquisition parameters, and taken by pointing to a number of
    sky positions, with separations of the order of sub arcsecs,
    either by moving the either by nodding the TS, tilting the TS
    M2 or shifting the EMIR DTU, the latter being the most likely
    option. Displacements are of the order of fraction of pixels.
    Images share the large majority of the sky positions so they can
    be coadded. Used for improving the spatial resolution of the
    resulting images and not valid for sky or superflat images.


    **Observing modes:**

        * Micro-dithered images


    Recipe to reduce observations obtained in imaging mode with microdithering.
    A critical piece of information
    here is a table that clearly specifies which images can be labelled as
    *science*, and which ones as *sky*. Note that some images are used both as
    *science* and *sky* (when the size of the targets are small compared to the
    offsets).

    **Observing modes:**
     * Micro-dithered images

    **Inputs:**


        * Offsets between them
        * Master Dark
        * Bad pixel mask (BPM)
        * Non-linearity correction polynomials
        * Master flat (twilight/dome flats)
        * Master background (thermal background, only in K band)
        * Detector model (gain, RN, lecture mode)
        * Average extinction in the filter
        * Astrometric calibration (TBD)

    **Outputs:**

     * Image with three extensions: final image scaled to the individual
       exposure time, variance  and exposure time map OR number of images
       combined (TBD).

    **Procedure:**

    Images are regridded to a integer subdivision of the pixel and then they
    are corrected from dark, non-linearity and flat. It should be desirable
    that the microdithering follows a pattern that can be easily translated to
    a subdivision of the pixel size (by an integer *n* = 2, 3, 4,...) that
    does not require a too high *n* value. An iterative process starts:

     * Sky is computed from each frame, using the list of sky images of each
       science frame. The objects are avoided using a mask (from the second
       iteration on).

     * The relatiev offsets are the nominal from the telescope. From the second
       iteration on, we refine them using bright objects.

     * We combine the sky-subtracted images, output is: a new image, a variance
       image and a exposure map/number of images used map.

     * An object mask is generated.

     * We recompute the sky map, using the object mask as an additional input.
       From here we iterate (typically 4 times).

     * Finally, the images are corrected from atmospheric extinction and flux
       calibrated.

     * A preliminary astrometric calibration can always be used (using
       the central coordinates of the pointing and the plate
       scale in the detector).
       A better calibration might be computed using available stars (TBD).

    """

    def run(self, ri):

        subpix = ri.subpixelization
        window = ri.window

        frame, catalog = self.process(ri,
                                      window=window, subpix=subpix,
                                      target_is_sky=True,
                                      stop_after=DirectImageCommon.FULLRED)

        result = self.create_result(frame=frame, catalog=catalog)
        return result
