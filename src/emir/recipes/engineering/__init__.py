
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

from numina.core import BaseRecipe, Parameter, provides, requires
from numina.core import DataFrame

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

#@provides(DTU_XY_Calibration)
class DTU_XY_CalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * DTU X_Y calibration
    '''

#    __requires__ = [
#        Parameter('slit_pattern', None, 'Slit pattern'),
#        Parameter('dtu_range', None, 'DTU range: begin, end and step')
#    ]

    def __init__(self):
        super(DTU_XY_CalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_XY_Calibration()]}
    
#@provides(DTU_Z_Calibration)
class DTU_Z_CalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * DTU Z calibration
    '''

#    __requires__ = [
#        Parameter('dtu_range', None, 'DTU range: begin, end and step')
#    ]

    def __init__(self):
        super(DTU_Z_CalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_Z_Calibration()]}

#@requires()
#@provides(DTUFlexureCalibration)
class DTUFlexureRecipe(BaseRecipe):
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

#@provides(DTU_XY_Calibration)
class CSU2DetectorRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * CSU2Detector calibration
    '''

#    __requires__ = [
#        Parameter('dtu_range', None, 'DTU range: begin, end and step')
#    ]

    def __init__(self):
        super(CSU2DetectorRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [DTU_XY_Calibration()]}

#@requires()
#@provides(PointingOriginCalibration)
class FocalPlaneCalibrationRecipe(BaseRecipe):
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

#@requires()
#@provides(WavelengthCalibration)
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

#@requires()
#@provides(PointingOriginCalibration)
class RotationCenterRecipe(BaseRecipe):
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

#@requires()
#@provides(DataFrame)
class AstrometricCalibrationRecipe(BaseRecipe):
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
        return {'products': [DataFrame(None)]}

#@provides(PhotometricCalibration)
class PhotometricCalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * Photometric calibration
    '''

#    __requires__ = [
#        Parameter('phot', None, 'Information about standard stars')
#    ]

    def __init__(self):
        super(PhotometricCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [PhotometricCalibration()]}

#@provides(SpectroPhotometricCalibration)
class SpectroPhotometricCalibrationRecipe(BaseRecipe):
    '''

    **Observing modes:**

        * Spectrophotometric calibration
    '''

#    __requires__ = [
#        Parameter('sphot', None, 'Information about standard stars')
#    ]

    def __init__(self):
        super(SpectroPhotometricCalibrationRecipe, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [SpectroPhotometricCalibration()]}



