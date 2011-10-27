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

import StringIO

from pyfits import Header

from numina.compatibility import get_data
from numina.image.storage import FITSCreator


_result_types = ['image', 'spectrum']
_extensions = ['primary', 'variance', 'map', 'wcs']


_table = {('image','primary'): 'image_primary.txt',
          ('image', 'map'): 'image_map.txt',
          ('image', 'wcs'): 'image_wcs.txt',
          ('image', 'variance'): 'image_variance.txt',
          ('spectrum','primary'): 'spectrum_primary.txt',
          ('spectrum', 'map'): 'image_map.txt',
          ('spectrum', 'wcs'): 'image_wcs.txt',
          ('spectrum', 'variance'): 'image_variance.txt',
          }

def load_header(res, ext):
    try:
        res = _table[(res, ext)]
    except KeyError:
        return Header()    
    sfile = StringIO.StringIO(get_data('emir.instrument', res))
    hh = Header(txtfile=sfile)
    return hh

def load_all_headers():
    result = {}
    for res in _result_types:
        result[res] = {}
        for ext in _extensions:
            result[res][ext] = load_header(res, ext)
    
    return result
    
default = load_all_headers()

class EmirImageCreator(FITSCreator):
    '''Builder of Emir direct image.'''
    def __init__(self): 
        super(EmirImageCreator, self).__init__(default['image'])
        
class EmirSpectrumCreator(FITSCreator):
    '''Builder of Emir spectrum.'''
    def __init__(self): 
        super(EmirSpectrumCreator, self).__init__(default['spectrum'])
