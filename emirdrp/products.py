#
# Copyright 2008-2018 Universidad Complutense de Madrid
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
import copy

import six
import numpy

from numina.ext.gtc import DF
from numina.core import DataFrameType, DataProductType
import numina.types.array as arrtype
from numina.types.linescatalog import LinesCatalog
import numina.exceptions
import numina.types.product as prodtypes
import numina.types.structured
import numina.types.obsresult as obtypes

import emirdrp.datamodel

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


class EmirFrame(DataFrameType):

    def __init__(self, *args, **kwargs):
        super(EmirFrame, self).__init__(datamodel=emirdrp.datamodel.EmirDataModel)


class ProcessedFrame(EmirFrame):
    """A processed frame"""
    pass


class ProcessedImage(ProcessedFrame):
    """A processed image"""
    pass


class ProcessedMOS(ProcessedFrame):
    """A processed image with slit spectra"""
    pass


class ProcessedSpectrum(ProcessedFrame):
    """A 1d spectrum"""
    pass


class ProcessedImageProduct(prodtypes.DataProductMixin, ProcessedImage):
    pass


class ProcessedMOSProduct(prodtypes.DataProductMixin, ProcessedMOS):
    pass


class ProcessedSpectrumProduct(prodtypes.DataProductMixin, ProcessedSpectrum):
    pass


class MasterBadPixelMask(ProcessedImageProduct):
    pass


class MasterBias(ProcessedImageProduct):
    """Master bias product
    """
    pass


class MasterDark(ProcessedImageProduct):
    """Master dark product
    """
    pass


class SkyLinesCatalog(LinesCatalog):
    """Sky Lines Catalog
    """
    pass


class MasterIntensityFlat(ProcessedImageProduct):
    __tags__ = ['filter']

class MasterSpectralFlat(ProcessedImageProduct):
    __tags__ = ['grism', 'filter']


# FIXME: This is not really a calibration
class MasterSky(ProcessedImageProduct):
    __tags__ = ['filter']


# FIXME: This is not really a calibration
class SkySpectrum(ProcessedImageProduct):
    __tags__ = ['grism', 'filter']


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


class CoordinateListNType(arrtype.ArrayNType):
    def __init__(self, dimensions, default=None):
        super(CoordinateListNType,
              self).__init__(dimensions,
                             default=default)

    def validate(self, obj):
        ndims = len(obj.shape)
        if ndims != 2:
            raise numina.exceptions.ValidationError(
                '%r is not a valid %r' %
                (obj, self.__class__.__name__)
            )
        if obj.shape[1] != self.N:
            raise numina.exceptions.ValidationError(
                '%r is not a valid %r' %
                (obj, self.__class__.__name__)
            )


class CoordinateList1DType(CoordinateListNType):
    def __init__(self, default=None):
        super(CoordinateList1DType, self).__init__(1, default=default)


class CoordinateList2DType(CoordinateListNType):
    def __init__(self, default=None):
        super(CoordinateList2DType, self).__init__(2, default=default)
        self.add_dialect_info('gtc', DF.TYPE_DOUBLE_ARRAY2D)


class NominalPositions(prodtypes.DataProductMixin, CoordinateList2DType):
    pass


def default_nominal_positions():
    """Read default value of bnp"""
    import pkgutil
    try:
        import StringIO as S
    except ImportError:
        import io as S
    bardata = pkgutil.get_data('emirdrp.instrument.configs', 'bars_nominal_positions_test.txt')
    ss = S.StringIO(bardata.decode('utf8'))
    bars_nominal_positions = numpy.loadtxt(ss)
    return bars_nominal_positions


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


class RefinedBoundaryModelParam(numina.types.structured.BaseStructuredCalibration):
    """Refined parameters of MOS model
    """
    def __init__(self, instrument='unknown'):
        super(RefinedBoundaryModelParam, self).__init__(instrument)
        self.tags = {
            'grism': "unknown",
            'filter': "unknown"
        }
        self.contents = []

    def __getstate__(self):
        state = super(RefinedBoundaryModelParam, self).__getstate__()
        if six.PY2:
            state['contents'] = copy.copy(self.contents)
        else:
            state['contents'] = self.contents.copy()
        return state

    def __setstate__(self, state):
        super(RefinedBoundaryModelParam, self).__setstate__(state)
        if six.PY2:
            self.contents = copy.copy(state['contents'])
        else:
            self.contents = state['contents'].copy()

    def tag_names(self):
        return ['grism', 'filter']


class RectWaveCoeff(numina.types.structured.BaseStructuredCalibration):
    """Rectification and Wavelength Calibration Coefficients
    """
    def __init__(self, instrument='unknown'):
        super(RectWaveCoeff, self).__init__(instrument)
        self.tags = {
            'grism': "unknown",
            'filter': "unknown"
        }
        self.global_integer_offset_x_pix = 0
        self.global_integer_offset_y_pix = 0
        self.total_slitlets = 0
        self.missing_slitlets = []
        self.contents = []

    def __getstate__(self):
        state = super(RectWaveCoeff, self).__getstate__()

        keys = ['global_integer_offset_x_pix', 'global_integer_offset_y_pix',
                'total_slitlets', 'missing_slitlets']
        for key in keys:
            state[key] = self.__dict__[key]

        if six.PY2:
            state['contents'] = copy.copy(self.contents)
        else:
            state['contents'] = self.contents.copy()
        return state

    def __setstate__(self, state):
        super(RectWaveCoeff, self).__setstate__(state)

        keys = ['global_integer_offset_x_pix', 'global_integer_offset_y_pix',
                'total_slitlets', 'missing_slitlets']
        for key in keys:
            self.__dict__[key] = state[key]
        if six.PY2:
            self.contents = copy.copy(state['contents'])
        else:
            self.contents = state['contents'].copy()

    def tag_names(self):
        return ['grism', 'filter']


class MasterRectWave(numina.types.structured.BaseStructuredCalibration):
    """Rectification and Wavelength Calibration Library Product
    """
    def __init__(self, instrument='unknown'):
        import numina.core.tagexpr as tagexpr

        datamodel = emirdrp.datamodel.EmirDataModel()
        super(MasterRectWave, self).__init__(instrument, datamodel=datamodel)
        self.tags = {
            'grism': "unknown",
            'filter': "unknown"
        }

        my_tag_table = self.datamodel.query_attrs
        objtags = [my_tag_table[t] for t in self.tag_names()]
        self.query_expr = tagexpr.query_expr_from_attr(objtags)
        self.names_t = self.query_expr.tags()
        self.names_f = self.query_expr.fields()
        self.query_opts = []

        self.total_slitlets = 0
        self.missing_slitlets = []
        self.contents = []

    def __getstate__(self):
        state = super(MasterRectWave, self).__getstate__()

        keys = ['total_slitlets', 'missing_slitlets']
        for key in keys:
            state[key] = self.__dict__[key]

        if six.PY2:
            # list has no .copy method in PY2
            state['contents'] = copy.copy(self.contents)
        else:
            state['contents'] = self.contents.copy()
        return state

    def __setstate__(self, state):
        from numina.types.typedialect import dialect_info
        import numina.core.tagexpr as tagexpr

        super(MasterRectWave, self).__setstate__(state)

        # Additional values
        # These values are not set with setstate
        self.datamodel = emirdrp.datamodel.EmirDataModel()
        self.internal_default = None
        self.internal_type = self.__class__
        self.internal_dialect = dialect_info(self)
        #
        my_tag_table = self.datamodel.query_attrs

        objtags = [my_tag_table[t] for t in self.tag_names()]
        self.query_expr = tagexpr.query_expr_from_attr(objtags)
        self.names_t = self.query_expr.tags()
        self.names_f = self.query_expr.fields()
        self.query_opts = []

        keys = ['total_slitlets', 'missing_slitlets']
        for key in keys:
            self.__dict__[key] = state[key]

        if six.PY2:
            # list has no .copy method in PY2
            self.contents = copy.copy(state['contents'])
        else:
            self.contents = state['contents'].copy()

    def tag_names(self):
        return ['grism', 'filter']

try:
    # FIXME: put this where it makes sense
    from gtc.DSL.DGCSTypes.IDL_Adapters import toIDL_, toIDL_struct, toElementType

    toIDL_.register(MasterRectWave, toIDL_struct)

    @toElementType.register(MasterRectWave)
    def _(value):
        return DF.TYPE_STRUCT
except ImportError:
    pass


if __name__ == '__main__':
    m = MasterRectWave().query_expr
    print(m.tags())
    m = MasterIntensityFlat().query_expr
    print(m.tags())
