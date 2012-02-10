
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

from numina.recipes import RecipeBase, Parameter, provides, requires
from numina.recipes import Image

from ...dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from ...dataproducts import TelescopeOffset, MSMPositions
from ...dataproducts import MasterIntensityFlat

from ...dataproducts import DTU_XY_Calibration, DTU_Z_Calibration
from ...dataproducts import CSU2DetectorCalibration, DTUFlexureCalibration
from ...dataproducts import (PointingOriginCalibration, 
                             SpectroPhotometricCalibration)
from ...dataproducts import PhotometricCalibration, WavelengthCalibration
__all__ = []

_logger = logging.getLogger('emir.recipes')

@provides(DTU_XY_Calibration)
class DTU_XY_CalibrationRecipe(RecipeBase):
    '''

    **Observing modes:**

        * DTU X_Y calibration
    '''

    __requires__ = [
        Parameter('slit_pattern', None, 'Slit pattern'),
        Parameter('dtu_range', None, 'DTU range: begin, end and step')
    ]

    def __init__(self):
        super(DTU_XY_CalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_XY_Calibration()]}
    
@provides(DTU_Z_Calibration)
class DTU_Z_CalibrationRecipe(RecipeBase):
    '''

    **Observing modes:**

        * DTU Z calibration
    '''

    __requires__ = [
        Parameter('dtu_range', None, 'DTU range: begin, end and step')
    ]

    def __init__(self):
        super(DTU_Z_CalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_Z_Calibration()]}

@requires()
@provides(DTUFlexureCalibration)
class DTUFlexureRecipe(RecipeBase):
    '''

    **Observing modes:**

        * DTU Flexure compensation
    '''

    # __requires__ = []

    def __init__(self):
        super(DTUFlexureRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTUFlexureCalibration()]}

@provides(DTU_XY_Calibration)
class CSU2DetectorRecipe(RecipeBase):
    '''

    **Observing modes:**

        * CSU2Detector calibration
    '''

    __requires__ = [
        Parameter('dtu_range', None, 'DTU range: begin, end and step')
    ]

    def __init__(self):
        super(CSU2DetectorRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_XY_Calibration()]}

@requires()
@provides(PointingOriginCalibration)
class FocalPlaneCalibrationRecipe(RecipeBase):
    '''

    **Observing modes:**

        * Lateral color
    '''

    # __requires__ = []

    def __init__(self):
        super(FocalPlaneCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [PointingOriginCalibration()]}

@requires()
@provides(WavelengthCalibration)
class SpectralCharacterizationRecipe(RecipeBase):
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

@requires()
@provides(PointingOriginCalibration)
class RotationCenterRecipe(RecipeBase):
    '''

    **Observing modes:**

        * Centre of rotation
    '''

    # __requires__ = []

    def __init__(self):
        super(RotationCenterRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [PointingOriginCalibration()]}

@requires()
@provides(Image)
class AstrometricCalibrationRecipe(RecipeBase):
    '''

    **Observing modes:**

        * Astrometric calibration
    '''

    # __requires__ = []

    def __init__(self):
        super(AstrometricCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [Image(None)]}

@provides(PhotometricCalibration)
class PhotometricCalibrationRecipe(RecipeBase):
    '''

    **Observing modes:**

        * Photometric calibration
    '''

    __requires__ = [
        Parameter('phot', None, 'Information about standard stars')
    ]

    def __init__(self):
        super(PhotometricCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [PhotometricCalibration()]}

@provides(SpectroPhotometricCalibration)
class SpectroPhotometricCalibrationRecipe(RecipeBase):
    '''

    **Observing modes:**

        * Spectrophotometric calibration
    '''

    __requires__ = [
        Parameter('sphot', None, 'Information about standard stars')
    ]

    def __init__(self):
        super(SpectroPhotometricCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [SpectroPhotometricCalibration()]}



