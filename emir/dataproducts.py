#
# Copyright 2008-2011 Sergio Pascual
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

from emir.instrument.headers import default

_logger = logging.getLogger('emir.dataproducts')

# Image -> Raw: PRIMARY
#       -> Result: PRIMARY, VARIANCE, MAP 

_result_types = ['image', 'spectrum']
_extensions = ['primary', 'variance', 'map', 'wcs']

def create_raw(data=None, headers=None, default_headers=None, imgtype='image'):
    
    if imgtype not in _result_types:
        raise TypeError('%s not in %s' % (imgtype, _result_types))
    
    hdefault = default_headers or default
    
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
    
    hdefault = default_headers or default
    
    for extname in ['variance', 'map']:
        edata = extensions[extname]
        hdu = pyfits.ImageHDU(edata[0], hdefault[imgtype][extname], name=extname)
        hdulist.append(hdu)
    return hdulist

