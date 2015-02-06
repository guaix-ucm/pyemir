
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

''' Recipes for engineering modes.

'''

import logging

from numina.core import Parameter
from numina.core import DataFrameType, Product, RecipeRequirements
from numina.core.requirements import ObservationResultRequirement

from emir.core import EmirRecipe
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import TelescopeOffset, MSMPositions
from emir.dataproducts import MasterIntensityFlat

from emir.dataproducts import DTU_XY_Calibration, DTU_Z_Calibration
from emir.dataproducts import CSU2DetectorCalibration, DTUFlexureCalibration
from emir.dataproducts import (PointingOriginCalibration,
                               SpectroPhotometricCalibration)
from emir.dataproducts import PhotometricCalibration, WavelengthCalibration

_logger = logging.getLogger('numina.recipes.emir')

from .detectorgain import GainRecipe1
from .cosmetics import CosmeticsRecipe


class DTU_XY_CalibrationRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * DTU X_Y calibration
    '''
    
    obresult = ObservationResultRequirement()
    slit_pattern = Parameter([], 'Slit pattern'),
    dtu_range = Parameter([], 'DTU range: begin, end and step')

    calibration = Product(DTU_XY_Calibration)


    def run(self, rinput):
        return self.create_result(calibration=DTU_XY_Calibration())


class DTU_Z_CalibrationRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * DTU Z calibration
    '''
    obresult = ObservationResultRequirement()
    dtu_range = Parameter([], 'DTU range: begin, end and step')

    calibration = Product(DTU_Z_Calibration)

    def run(self, rinput):
        return self.create_result(calibration=DTU_Z_Calibration())


class DTUFlexureRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * DTU Flexure compensation
    '''

    obresult = ObservationResultRequirement()
    calibration = Product(DTUFlexureCalibration)

    def run(self, rinput):
        return self.create_result(calibration=DTUFlexureCalibration())


class CSU2DetectorRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * CSU2Detector calibration
    '''

    obresult = ObservationResultRequirement()
    dtu_range = Parameter([], 'DTU range: begin, end and step')

    calibration = Product(DTU_XY_Calibration)

    def run(self, rinput):
        return self.create_result(calibration=DTU_XY_Calibration())


class FocalPlaneCalibrationRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * Lateral color
    '''

    obresult = ObservationResultRequirement()

    calibration = Product(PointingOriginCalibration)


    def run(self, rinput):
        return self.create_result(calibration=PointingOriginCalibration())


class SpectralCharacterizationRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * Spectral characterization
    '''

    obresult = ObservationResultRequirement()
    calibration = Product(WavelengthCalibration)

    def run(self, rinput):
        return self.create_result(calibration=WavelengthCalibration())


class RotationCenterRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * Centre of rotation
    '''

    obresult = ObservationResultRequirement()
    calibration = Product(PointingOriginCalibration)

    def run(self, rinput):
        return self.create_result(calibration=PointingOriginCalibration())


class AstrometricCalibrationRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * Astrometric calibration
    '''

    obresult = ObservationResultRequirement()
    calibration = Product(DataFrameType)

    def run(self, rinput):
        return self.create_result(calibration=None)


class PhotometricCalibrationRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * Photometric calibration
    '''

    obresult = ObservationResultRequirement()
    phot = Parameter([], 'Information about standard stars')
    calibration = Product(PhotometricCalibration)


    def run(self, rinput):
        return self.create_result(calibration=PhotometricCalibration())


class SpectroPhotometricCalibrationRecipe(EmirRecipe):

    '''

    **Observing modes:**

        * Spectrophotometric calibration
    '''

    obresult = ObservationResultRequirement()
    sphot = Parameter([], 'Information about standard stars')

    calibration = Product(SpectroPhotometricCalibration)

    def run(self, rinput):
        return self.create_result(calibration=SpectroPhotometricCalibration())
