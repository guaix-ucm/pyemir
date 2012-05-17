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

import logging
import pyfits

from numina.recipes import DataFrame, DataProduct

from .simulator import EmirImageFactory

_logger = logging.getLogger('emir.dataproducts')

class MasterBadPixelMask(DataFrame):
    def __init__(self, hdu):
        super(MasterBias, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['ccdmode']

class MasterBias(DataFrame):
    '''Master bias product
    
    This image has 4 extensions: primary, two variance extensions
    and number of pixels used in the combination.
    
    The variance extensions are computed using two different methods. 
    The first one is the variance of the same pixels in different images.
    The second extension is the variance of each channel in the final image.
    
    
    '''
    def __init__(self, hdu):
        super(MasterBias, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['ccdmode']

class MasterDark(DataFrame):
    '''Master dark product
    
    This image has 4 extensions: primary, two variance extensions
    and number of pixels used in the combination.
    
    The variance extensions are computed using two different methods. 
    The first one is the variance of the same pixels in different images.
    The second extension is the variance of each channel in the final image.
    
    
    '''
    def __init__(self, hdu):
        super(MasterDark, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['ccdmode']

class MasterIntensityFlat(DataFrame):
    def __init__(self, hdu):
        super(MasterIntensityFlat, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['ccdmode']
        yield 'filter0', hdr['filter']
        
class MasterSpectralFlat(DataFrame):
    def __init__(self, hdu):
        super(MasterSpectralFlat, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['ccdmode']
        yield 'filter0', hdr['filter']

class Spectra(DataFrame):
    def __init__(self, hdu):
        super(Spectra, self).__init__(hdu)

    def metadata(self):
        hdr = self.image[0].header
        yield 'detector0.mode', hdr['ccdmode']
        yield 'filter0', hdr['filter']
        
class DataCube(DataFrame):
    def __init__(self, hdu):
        super(DataFrame, self).__init__(None)

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
    pass

class MasterRONMap(DataProduct):
    pass

class NonLinearityCalibration(DataProduct):
    def __init__(self, poly):
        super(NonLinearityCalibration, self).__init__()
        self.poly = poly

class TelescopeOffset(DataProduct):
    pass

class MSMPositions(DataProduct):
    pass

class SourcesCatalog(DataProduct):
    pass

class LinesCatalog(DataProduct):
    pass

class ChannelLevelStatistics(DataProduct):
    ''' A list of exposure time, mean, std dev and median per channel'''
    def __init__(self, exposure):
        self.exposure = exposure
        self.statistics = []

# DataFrame -> Raw: PRIMARY
#       -> Result: PRIMARY, VARIANCE, MAP 

_result_types = ['image', 'spectrum']
_extensions = ['primary', 'variance', 'map', 'wcs']

def create_raw(data=None, headers=None, default_headers=None, imgtype='image'):
    
    if imgtype not in _result_types:
        raise TypeError('%s not in %s' % (imgtype, _result_types))
    
    hdefault = default_headers or EmirImageFactory.default
    
    hdu = pyfits.PrimaryHDU(data, hdefault[imgtype]['primary'])
                
    if headers is not None:
        _logger.info('Updating keywords in %s header', 'PRIMARY')      
        for key in headers:
            _logger.debug('Updating keyword %s with value %s', 
                          key, headers[key])
            hdu.header.update(key, headers[key])
    return pyfits.HDUList([hdu])


def create_result(data=None, headers=None, 
                   variance=None, variance_headers=None,
                   exmap=None, exmap_headers=None,
                   default_headers=None, imgtype='image'):
    hdulist = create_raw(data, headers, default_headers, imgtype)
    extensions = {}
    extensions['variance'] = (variance, variance_headers)
    extensions['map'] = (exmap, exmap_headers)
    
    hdefault = default_headers or EmirImageFactory.default
    
    for extname in ['variance', 'map']:
        edata = extensions[extname]
        hdu = pyfits.ImageHDU(edata[0], hdefault[imgtype][extname], name=extname)
        hdulist.append(hdu)
    return hdulist

