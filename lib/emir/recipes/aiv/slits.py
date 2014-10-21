#
# Copyright 2013-2014 Universidad Complutense de Madrid
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

'''Recipe to detect slits in the AIV mask'''

from __future__ import division
#
import logging
# import math
#
# import numpy
# import scipy.interpolate as itpl
# import scipy.optimize as opz
# from astropy.modeling import models, fitting
from astropy.io import fits
# import photutils
#
# from numina.array.recenter import img_box, centering_centroid
#
from numina.core import BaseRecipe, RecipeRequirements, RecipeError
from numina.core import Requirement, Product, DataProductRequirement, Parameter
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement
#
from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import DataFrameType, MasterIntensityFlat
# from emir.dataproducts import CoordinateList2DType
# from emir.dataproducts import ArrayType

# from photutils import aperture_circular
#
# from .procedures import compute_fwhm_spline2d_fit
# from .procedures import compute_fwhm_enclosed_direct
# from .procedures import compute_fwhm_enclosed_grow
# from .procedures import compute_fwhm_simple
# from .procedures import moments
# from .procedures import AnnulusBackgroundEstimator
# from .procedures import img_box2d
from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs 

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"

GAUSS_FWHM_FACTOR = 2.354800


class TestSlitDetectionRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')
    master_sky = DataProductRequirement(MasterIntensityFlat, 'Master Sky calibration')
#    pinhole_nominal_positions = Requirement(CoordinateList2DType, 'Nominal positions of the pinholes')
#    shift_coordinates = Parameter(True, 'Use header information to shift the pinhole positions from (0,0) to X_DTU, Y_DTU')
#    box_half_size = Parameter(4, 'Half of the computation box size in pixels')
#    recenter = Parameter(True, 'Recenter the pinhole coordinates')
#    max_recenter_radius = Parameter(2.0, 'Maximum distance for recentering')


class TestSlitDetectionRecipeResult(RecipeResult):
    frame = Product(DataFrameType)
#    positions = Product(ArrayType)
#    positions_alt = Product(ArrayType)
#    DTU = Product(ArrayType)
#    filter = Product(str)
#    readmode = Product(str)
#    IPA = Product(float)


@define_requirements(TestSlitDetectionRecipeRequirements)
@define_result(TestSlitDetectionRecipeResult)
class TestSlitDetectionRecipe(BaseRecipe):

    def __init__(self):
        super(TestSlitDetectionRecipe, self).__init__(
            author=_s_author,
            version="0.1.0"
            )

    def run(self, rinput):
        _logger.info('starting pinhole processing')


        flow = init_filters_bdfs(rinput)
        
        hdu = basic_processing_with_combination(rinput, flow=flow)
        
        hdr = hdu.header
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')

        _logger.debug('finding pinholes')

        try:
            filtername = hdr['FILTER']
            readmode = hdr['READMODE']
            ipa = hdr['IPA']
            xdtu = hdr['XDTU']
            ydtu = hdr['YDTU']
            zdtu = hdr['ZDTU']
        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

        result = self.create_result(frame=hdu)

        return result
