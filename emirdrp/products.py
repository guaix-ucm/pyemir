#
# Copyright 2008-2017 Universidad Complutense de Madrid
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

"""Data products produced by the EMIR pipeline."""

import uuid

import numpy

from numina.ext.gtc import DF
from numina.core import DataFrameType, DataProductType
from numina.core.products import ArrayNType
from numina.core.products import DataProductTag
from numina.core.requirements import InstrumentConfigurationType
from numina.core import ValidationError


base_schema_description = {
    'keywords': {
        'INSTRUME': {'mandatory': True},
        'READMODE': {'mandatory': True,
                     'value': ['SIMPLE', 'BIAS', 'SINGLE',
                               'CDS', 'FOWLER', 'RAMP']
                     },
        'EXPTIME': {'value': float},
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
                     'value': ['SIMPLE', 'BIAS', 'SINGLE',
                               'CDS', 'FOWLER', 'RAMP']
                     },
        'BUNIT': {'value': ['ADU', 'ADU/s']},
        'IMGTYPE': {'mandatory': True, 'type': 'string'},
        }
    }


class EMIRImageProduct(DataProductTag, DataFrameType):

    def convert_out(self, obj):
        newobj = super(EMIRImageProduct, self).convert_out(obj)
        if newobj:
            hdulist = newobj.open()
            hdr = hdulist[0].header
            if 'EMIRUUID' not in hdr:
                hdr['EMIRUUID'] = str(uuid.uuid1())
        return newobj


class EMIRConfigurationType(InstrumentConfigurationType):

    def validate(self, value):
        super(EMIRConfigurationType, self).validate(value)


class MasterBadPixelMask(EMIRImageProduct):
    pass


class MasterBias(EMIRImageProduct):
    """Master bias product

    This image has 4 extensions: primary, two variance extensions
    and number of pixels used in the combination.

    The variance extensions are computed using two different methods.
    The first one is the variance of the same pixels in different images.
    The second extension is the variance of each channel in the final image.
    """
    pass


class MasterDark(EMIRImageProduct):
    """Master dark product

    This image has 4 extensions: primary, two variance extensions
    and number of pixels used in the combination.

    The variance extensions are computed using two different methods.
    The first one is the variance of the same pixels in different images.
    The second extension is the variance of each channel in the final image.
    """
    pass



class MasterIntensityFlat(EMIRImageProduct):
    pass


class MasterSky(EMIRImageProduct):
    pass


class SkySpectrum(EMIRImageProduct):
    pass


class MasterSpectralFlat(EMIRImageProduct):
    pass


class Spectra(DataFrameType):
    pass


class DataCube(DataFrameType):
    pass


class TelescopeFocus(DataProductType):
    def __init__(self, default=None):
        super(TelescopeFocus, self).__init__(ptype=float)


class DTUFocus(DataProductType):
    def __init__(self, default=None):
        super(DTUFocus, self).__init__(ptype=float)


class DTU_XY_Calibration(DataFrameType):
    pass


class DTU_Z_Calibration(DataFrameType):
    pass


class DTUFlexureCalibration(DataFrameType):
    pass


# FIXME:
class SlitTransmissionCalibration(DataFrameType):
    pass


# FIXME:
class WavelengthCalibration(DataFrameType):
    pass


class CSU2DetectorCalibration(DataFrameType):
    pass


class PointingOriginCalibration(DataFrameType):
    pass


class SpectroPhotometricCalibration(DataFrameType):
    pass


class PhotometricCalibration(DataFrameType):
    pass


class MasterGainMap(DataProductType):
    def __init__(self, mean, var, frame):
        self.mean = mean
        self.var = var
        self.frame = frame

    def __getstate__(self):
        gmean = list(map(float, self.mean.flat))
        gvar = list(map(float, self.var.flat))
        return {'frame': self.frame, 'mean': gmean, 'var': gvar}


class MasterRONMap(DataProductType):
    def __init__(self, mean, var):
        self.mean = mean
        self.var = var

    def __getstate__(self):
        gmean = map(float, self.mean.flat)
        gvar = map(float, self.var.flat)
        return {'mean': gmean, 'var': gvar}


class TelescopeOffset(DataProductType):
    def __init__(self, default=None):
        super(TelescopeOffset, self).__init__(ptype=float)


class CoordinateListNType(ArrayNType):
    def __init__(self, dimensions, default=None):
        super(CoordinateListNType,
              self).__init__(dimensions,
                             default=default)

    def validate(self, obj):
        ndims = len(obj.shape)
        if ndims != 2:
            raise ValidationError('%r is not a valid %r' %
                                  (obj, self.__class__.__name__)
                                  )
        if obj.shape[1] != self.N:
            raise ValidationError('%r is not a valid %r' %
                                  (obj, self.__class__.__name__)
                                  )


class CoordinateList1DType(CoordinateListNType):
    def __init__(self, default=None):
        super(CoordinateList1DType, self).__init__(1, default=default)


class CoordinateList2DType(CoordinateListNType):
    def __init__(self, default=None):
        super(CoordinateList2DType, self).__init__(2, default=default)
        self.add_dialect_info('gtc', DF.TYPE_DOUBLE_ARRAY2D)


class NominalPositions(DataProductTag, CoordinateList2DType):
    pass


class MSMPositions(DataProductType):
    def __init__(self):
        super(MSMPositions, self).__init__(ptype=list)


class SourcesCatalog(DataProductType):
    def __init__(self):
        super(SourcesCatalog, self).__init__(ptype=list)


class CentroidsTableType(DataProductType):
    """Table with information about focus centroids."""
    def __init__(self):
        super(CentroidsTableType, self).__init__(ptype=numpy.ndarray)


class ChannelLevelStatistics(DataProductType):
    """A list of exposure time, mean, std dev and median per channel"""
    def __init__(self, exposure, statistics):
        self.exposure = exposure
        self.statistics = statistics

    def __numina_dump__(self, obj, where):
        fname = 'statistics.txt'

        header = ("Channel Level Statistics\n"
                  "comment 2\n"
                  "pixels start in 1\n"
                  "pixels end in 2048\n"
                  "exposure={exposure}\n"
                  "xbegin xend ybegin yend mean median var\n"
                  )

        inter = header.format(exposure=obj.exposure)
        numpy.savetxt(fname, obj.statistics, header=inter)
        return fname



class ChannelLevelStatisticsType(DataProductType):
    """A list of exposure time, mean, std dev and median per channel"""
    def __init__(self):
        super(ChannelLevelStatisticsType,
              self).__init__(ptype=ChannelLevelStatistics)


class LinesCatalog(DataProductType):
    def __init__(self):
        super(LinesCatalog, self).__init__(ptype=numpy.ndarray)
