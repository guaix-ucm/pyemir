
# Copyright 2008-2012 Universidad Complutense de Madrid
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

from numina.core import BaseRecipe, Parameter, define_result, define_input
from numina.core import FrameDataProduct, RecipeResult, Product, RecipeInput 

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

class DTU_XY_CalibrationRecipeInput(RecipeInput):
    slit_pattern = Parameter(None, 'Slit pattern'),
    dtu_range = Parameter(None, 'DTU range: begin, end and step')

class DTU_XY_CalibrationRecipeResult(RecipeResult):
    calibration = Product(DTU_XY_Calibration)

@define_input(DTU_XY_CalibrationRecipeInput)
@define_result(DTU_XY_CalibrationRecipeResult)
class DTU_XY_CalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * DTU X_Y calibration
    '''

    def __init__(self):
        super(DTU_XY_CalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_XY_Calibration()]}

class DTU_Z_CalibrationRecipeInput(RecipeInput):
    dtu_range = Parameter(None, 'DTU range: begin, end and step')

class DTU_Z_CalibrationRecipeResult(RecipeResult):
    calibration = Product(DTU_Z_Calibration)

@define_input(DTU_Z_CalibrationRecipeInput)
@define_result(DTU_Z_CalibrationRecipeResult) 
class DTU_Z_CalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * DTU Z calibration
    '''

    def __init__(self):
        super(DTU_Z_CalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_Z_Calibration()]}

class DTU_FlexureRecipeInput(RecipeInput):
    pass

class DTU_FlexureRecipeResult(RecipeResult):
    calibration = Product(DTUFlexureCalibration)

@define_input(DTU_FlexureRecipeInput)
@define_result(DTU_FlexureRecipeResult)
class DTUFlexureRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * DTU Flexure compensation
    '''

    def __init__(self):
        super(DTUFlexureRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTUFlexureCalibration()]}

class CSU2DetectorRecipeInput(RecipeInput):
    dtu_range = Parameter(None, 'DTU range: begin, end and step')

class CSU2DetectorRecipeResult(RecipeResult):
    calibration = Product(DTU_XY_Calibration)

@define_input(CSU2DetectorRecipeInput)
@define_result(CSU2DetectorRecipeResult)
class CSU2DetectorRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * CSU2Detector calibration
    '''

    def __init__(self):
        super(CSU2DetectorRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_XY_Calibration()]}

class FocalPlaneCalibrationRecipeInput(RecipeInput):
    pass

class FocalPlaneCalibrationRecipeResult(RecipeResult):
    calibration = Product(PointingOriginCalibration)

@define_input(FocalPlaneCalibrationRecipeInput)
@define_result(FocalPlaneCalibrationRecipeResult)
class FocalPlaneCalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * Lateral color
    '''

    def __init__(self):
        super(FocalPlaneCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [PointingOriginCalibration()]}

class SpectralCharacterizationRecipeInput(RecipeInput):
    pass

class SpectralCharacterizationRecipeResult(RecipeResult):
    calibration = Product(WavelengthCalibration)

@define_input(SpectralCharacterizationRecipeInput)
@define_result(SpectralCharacterizationRecipeResult)
class SpectralCharacterizationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * Spectral characterization
    '''

    def __init__(self):
        super(SpectralCharacterizationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [WavelengthCalibration()]}

class RotationCenterRecipeInput(RecipeInput):
    pass

class RotationCenterRecipeResult(RecipeResult):
    calibration = Product(PointingOriginCalibration)

@define_input(RotationCenterRecipeInput)
@define_result(RotationCenterRecipeResult)
class RotationCenterRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * Centre of rotation
    '''

    def __init__(self):
        super(RotationCenterRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [PointingOriginCalibration()]}

class AstrometricCalibrationRecipeInput(RecipeInput):
    pass

class AstrometricCalibrationRecipeResult(RecipeResult):
    calibration = Product(FrameDataProduct)

@define_input(AstrometricCalibrationRecipeInput)
@define_result(AstrometricCalibrationRecipeResult)
class AstrometricCalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * Astrometric calibration
    '''

    def __init__(self):
        super(AstrometricCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return AstrometricCalibrationRecipeResult(calibration=None)

class PhotometricCalibrationRecipeInput(RecipeInput):
    phot = Parameter(None, 'Information about standard stars')
    
class PhotometricCalibrationRecipeResult(RecipeResult):
    calibration = Product(PhotometricCalibration)
    
@define_input(PhotometricCalibrationRecipeInput)
@define_result(PhotometricCalibrationRecipeResult)
class PhotometricCalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * Photometric calibration
    '''

    def __init__(self):
        super(PhotometricCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [PhotometricCalibration()]}

class SpectroPhotometricCalibrationRecipeInput(RecipeInput):
    sphot = Parameter(None, 'Information about standard stars')
    
class SpectroPhotometricCalibrationRecipeResult(RecipeResult):
    calibration = Product(SpectroPhotometricCalibration)
    
@define_input(SpectroPhotometricCalibrationRecipeInput)
@define_result(SpectroPhotometricCalibrationRecipeResult)
class SpectroPhotometricCalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * Spectrophotometric calibration
    '''

    def __init__(self):
        super(SpectroPhotometricCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [SpectroPhotometricCalibration()]}
