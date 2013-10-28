#
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

'''Data products produced by the EMIR pipeline.'''

from numina.core import FrameDataProduct, DataProduct

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

class DarkCurrentValue(DataProduct):
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

class DTU_XY_Calibration(DataProduct):
    pass

class DTU_Z_Calibration(DataProduct):
    pass

class DTUFlexureCalibration(DataProduct):
    pass

class SlitTransmissionCalibration(DataProduct):
    pass

class WavelengthCalibration(DataProduct):
    pass

class CSU2DetectorCalibration(DataProduct):
    pass

class PointingOriginCalibration(DataProduct):
    pass

class SpectroPhotometricCalibration(DataProduct):
    pass

class PhotometricCalibration(DataProduct):
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

class NonLinearityCalibration(DataProduct):
    def __init__(self, poly):
        super(NonLinearityCalibration, self).__init__(default=poly)
        self.poly = poly

class TelescopeOffset(DataProduct):
    pass

class MSMPositions(DataProduct):
    pass

class SourcesCatalog(DataProduct):
    pass

class LinesCatalog(DataProduct):
    pass

def _s_to_f(myslice):
    b = myslice.start
    e = myslice.stop
    return slice(b+1, e)

class ChannelLevelStatistics(DataProduct):
    ''' A list of exposure time, mean, std dev and median per channel'''
    def __init__(self, exposure=None):
        self.exposure = exposure
        self.statistics = []

    def __getstate__(self):
        fname = 'statistics.txt'
        fs = '{1.start:5} {1.stop:5} {0.start:5} {0.stop:5} {2[0]} {2[1]} {2[2]}\n'
        with open(fname, 'w') as fd:
            fd.write('# Channel Level Statistics\n')
            fd.write('# comment 2\n')
            fd.write('# pixels start in 1\n')
            fd.write('# pixels end in 2048\n')
            fd.write('## exposure=%s\n' % self.exposure)
            fd.write('# xbegin xend ybegin yend mean median var\n')
            fd.write('#\n')
            for (regy, regx), vals in self.statistics:
                fd.write(fs.format(_s_to_f(regy), _s_to_f(regx), vals))
        return dict(filename=fname)

