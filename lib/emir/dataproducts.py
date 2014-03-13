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

'''Data products produced by the EMIR pipeline.'''

import numpy

from numina.core import FrameDataProduct, DataProduct
from numina.core.requirements import InstrumentConfigurationType
    
    
    
class EMIRConfigurationType(InstrumentConfigurationType):
    
    def validate(self, value):
        super(EMIRConfigurationType, self).validate(value)

class MasterBadPixelMask(FrameDataProduct):
    pass

class MasterBias(FrameDataProduct):
    '''Master bias product
    
    This image has 4 extensions: primary, two variance extensions
    and number of pixels used in the combination.
    
    The variance extensions are computed using two different methods. 
    The first one is the variance of the same pixels in different images.
    The second extension is the variance of each channel in the final image.
    
    
    '''
    pass

class MasterDark(FrameDataProduct):
    '''Master dark product
    
    This image has 4 extensions: primary, two variance extensions
    and number of pixels used in the combination.
    
    The variance extensions are computed using two different methods. 
    The first one is the variance of the same pixels in different images.
    The second extension is the variance of each channel in the final image.
    
    
    '''
    pass

class DarkCurrentValue(FrameDataProduct):
    pass

class MasterIntensityFlat(FrameDataProduct):
    pass
        
class MasterSpectralFlat(FrameDataProduct):
    pass

class Spectra(FrameDataProduct):
    pass
 
class DataCube(FrameDataProduct):
    pass

class TelescopeFocus(DataProduct):
    pass

class DTUFocus(DataProduct):
    pass

class DTU_XY_Calibration(FrameDataProduct):
    pass

class DTU_Z_Calibration(FrameDataProduct):
    pass

class DTUFlexureCalibration(FrameDataProduct):
    pass

# FIXME:
class SlitTransmissionCalibration(FrameDataProduct):
    pass

# FIXME:
class WavelengthCalibration(FrameDataProduct):
    pass

class CSU2DetectorCalibration(FrameDataProduct):
    pass

class PointingOriginCalibration(FrameDataProduct):
    pass

class SpectroPhotometricCalibration(FrameDataProduct):
    pass

class PhotometricCalibration(FrameDataProduct):
    pass

class MasterGainMap(DataProduct):
    def __init__(self, mean, var, frame):
        self.mean = mean
        self.var = var
        self.frame = frame

    def __getstate__(self):
        gmean = map(float, self.mean.flat)
        gvar = map(float, self.var.flat)
        return {'frame': self.frame, 'mean': gmean, 'var': gvar}

class MasterRONMap(DataProduct):
    def __init__(self, mean, var):
        self.mean = mean
        self.var = var

    def __getstate__(self):
        gmean = map(float, self.mean.flat)
        gvar = map(float, self.var.flat)
        return {'mean': gmean, 'var': gvar}

    pass

class CoordinateListType(DataProduct):
    def __init__(self, default=None):
        super(CoordinateListType, self).__init__(ptype=numpy.ndarray, default=default)

class TelescopeOffset(DataProduct):
    pass

class MSMPositions(DataProduct):
    pass

class SourcesCatalog(DataProduct):
    def __init__(self):
        super(SourcesCatalog, self).__init__(ptype=list)


class LinesCatalog(DataProduct):
    pass

class ChannelLevelStatistics(DataProduct):
    ''' A list of exposure time, mean, std dev and median per channel'''
    def __init__(self, exposure, statistics):
        self.exposure = exposure
        self.statistics = statistics

class ChannelLevelStatisticsType(DataProduct):
    ''' A list of exposure time, mean, std dev and median per channel'''
    def __init__(self):
        super(ChannelLevelStatisticsType, self).__init__(ptype=ChannelLevelStatistics)

