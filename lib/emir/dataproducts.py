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
from numina.core import ValidationError

# FIXME
try:
    import ext.gtc
except ImportError:
    # We are not in GTC
    pass

base_schema_description = {
    'keywords': {
        'INSTRUME': {'mandatory': True},
        'READMODE': {'mandatory': True, 
            'value': ['SIMPLE', 'BIAS', 'SINGLE', 'CDS', 'FOWLER', 'RAMP']},
        'EXPTIME': {'value': float},
 #       'XDTU': {'mandatory': True, 'value': float},
 #       'YDTU': {'mandatory': True, 'value': float},
 #       'ZDTU': {'mandatory': True, 'value': float},
        'NUMINAID': {'value': int}
        }
    }

gtc_proc_schema_description = {
    'keywords': {
        'NUMINAID': {'mandatory': True, 'value': int}
        }
    }

emir_schema_description = {
    'keywords': {
        'INSTRUME': {'mandatory': True, 'value': 'EMIR'},
        'READMODE': {'mandatory': True, 
            'value': ['SIMPLE', 'BIAS', 'SINGLE', 'CDS', 'FOWLER', 'RAMP']},
        'BUNIT': {'value': ['ADU', 'ADU/s']},
        'IMGTYPE': {'mandatory': True, 'type': 'string'},
 #       'XDTU': {'mandatory': True, 'value': float},
 #       'YDTU': {'mandatory': True, 'value': float},
 #       'ZDTU': {'mandatory': True, 'value': float}, 
        }
    }

class MasterFrameProduct(FrameDataProduct):
    
    def __init__(self):
        super(MasterFrameProduct, self).__init__()
        self.headerschema.extend(gtc_proc_schema_description)
        

class EMIRConfigurationType(InstrumentConfigurationType):
    
    def validate(self, value):
        super(EMIRConfigurationType, self).validate(value)

class EMIRFrame(FrameDataProduct):
    
    def __init__(self):
        super(EMIRFrame, self).__init__()
        self.headerschema.extend(base_schema_description)
        self.headerschema.extend(emir_schema_description)
            
    def validate_hdu(self, hdu):
        self.headerschema.validate(hdu.header)
        
    def validate_hdulist(self, hdulist):
        super(EMIRFrame, self).validate_hdulist(hdulist)
        self.validate_hdu(hdulist[0])


class MasterBadPixelMask(EMIRFrame):
    pass

class RawBias(EMIRFrame):
    '''Raw bias frame'''

    def validate_hdu(self, hdu):
        super(RawBias, self).validate_hdu(hdu)
        # Check READMODE is valid
        header = hdu.header
        if header['READMODE'] not in ['SIMPLE', 'BIAS', 'SINGLE']:
            raise ValidationError('not a bias')
        
        return True

class RawDark(EMIRFrame):
    '''Raw dark frame'''
    def validate_hdu(self, hdu):
        super(RawDark, self).validate_hdu(hdu)

        header = hdu.header
        if header['IMGTYPE'] != 'DARK':
            raise ValidationError('not a dark')
        
        return True

class RawIntensityFlat(EMIRFrame):
    def validate_hdu(self, hdu):
        super(RawIntensityFlat, self).validate_hdu(hdu)
        return True

class MasterBias(RawBias, MasterFrameProduct):
    '''Master bias product
    
    This image has 4 extensions: primary, two variance extensions
    and number of pixels used in the combination.
    
    The variance extensions are computed using two different methods. 
    The first one is the variance of the same pixels in different images.
    The second extension is the variance of each channel in the final image.
    '''
    pass


class MasterDark(RawDark, MasterFrameProduct):
    '''Master dark product
    
    This image has 4 extensions: primary, two variance extensions
    and number of pixels used in the combination.
    
    The variance extensions are computed using two different methods. 
    The first one is the variance of the same pixels in different images.
    The second extension is the variance of each channel in the final image.
    '''
    pass

class DarkCurrentValue(EMIRFrame):
    pass

class MasterIntensityFlat(RawIntensityFlat, MasterFrameProduct):
    pass
        
class MasterSpectralFlat(EMIRFrame):
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

class TelescopeOffset(DataProduct):
    pass


class ArrayType(DataProduct):
    def __init__(self, default=None):
        super(ArrayType, self).__init__(ptype=numpy.ndarray, default=default)

    def store(self, obj):
        result = numpy.array(obj)
        return result


class CoordinateListNType(DataProduct):
    def __init__(self, dimensions, default=None):
        super(CoordinateListNType, self).__init__(ptype=numpy.ndarray, default=default)
        self.N = dimensions

    def validate(self, obj):
        ndims = len(obj.shape)
        if ndims != 2:
            raise ValidationError('%r is not a valid %r' % (obj, self.__class__.__name__))
        if obj.shape[1] != self.N:
            raise ValidationError('%r is not a valid %r' % (obj, self.__class__.__name__))
            

    def store(self, obj):
        result = numpy.array(obj)
        return result

class CoordinateList2DType(CoordinateListNType):
    def __init__(self, default=None):
        super(CoordinateList2DType, self).__init__(2, default=default)

class MSMPositions(DataProduct):
    pass

class SourcesCatalog(DataProduct):
    def __init__(self):
        super(SourcesCatalog, self).__init__(ptype=list)


class LinesCatalog(DataProduct):
    pass

class CentroidsTableType(DataProduct):
    '''Table with information about focus centroids.'''
    def __init__(self):
        super(CentroidsTableType, self).__init__(ptype=numpy.ndarray)

class ChannelLevelStatistics(DataProduct):
    ''' A list of exposure time, mean, std dev and median per channel'''
    def __init__(self, exposure, statistics):
        self.exposure = exposure
        self.statistics = statistics

class ChannelLevelStatisticsType(DataProduct):
    ''' A list of exposure time, mean, std dev and median per channel'''
    def __init__(self):
        super(ChannelLevelStatisticsType, self).__init__(ptype=ChannelLevelStatistics)

