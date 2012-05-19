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

'''EMIR simulator.'''

from pkgutil import get_data
from StringIO import StringIO
import logging

import pyfits
import numpy
from numina.treedict import TreeDict
from numina.instrument import Shutter
from numina.instrument.template import interpolate

from emir.instrument.detector import EmirDetector, EmirDas

_logger = logging.getLogger('emir.simulator')

class BaseInstrument(object):
    def __init__(self, detector):

        self.shutter = Shutter()

        self.detector = detector
        self.das = EmirDas(detector)
        
        # Connect components
        self.detector.connect(self.shutter)
        
        self.detector.meta['gain'] = 2.8
        self.detector.meta['readnoise'] = 3.5

        self.meta = TreeDict()
        self.meta['name'] = 'EMIR'
        self.meta['focalstation'] = 'NASMYTH'
        self.meta['filter'] = 'J'
        self.meta['imagetype'] = ''
        self.meta['detector'] = self.detector.meta
        self.meta['das'] = self.das.meta
    
    def connect(self, source):
        self.shutter.connect(source)
 
    def filter(self, name):
        self.meta['filter'] = name

    def acquire(self, time):
        img = self.das.acquire(time)
        return self.meta, img

    def imagetype(self, name):
        self.meta['imagetype'] = name

class Instrument(BaseInstrument):
    def __init__(self):
        
        # FIXME: Testing flat
        xx, yy= numpy.mgrid[-1:1:2048j, -1:1:2048j]
        flat = 2 - 0.7 * (xx**2+yy**2)
        flat /= flat.mean()
        
        detector = EmirDetector(flat=flat)
        BaseInstrument.__init__(self, detector)

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

def _load_header(res, ext):
    try:
        res = _table[(res, ext)]
    except KeyError:
        return pyfits.Header()    
    sfile = StringIO(get_data('emir.instrument', res))
    hh = pyfits.Header(txtfile=sfile)
    return hh

def _load_all_headers():
    result = {}
    for res in _result_types:
        result[res] = {}
        for ext in _extensions:
            result[res][ext] = _load_header(res, ext)
    
    return result
    
class EmirImageFactory(object):
    
    default = _load_all_headers() 

    # FIXME: workaround
    def __init__(self):
        sfile = StringIO(get_data('emir.instrument', 'image_primary.txt'))
        self.p_templ = pyfits.Header(txtfile=sfile)
        del sfile

    def create(self, meta, data):
        hh = self.p_templ.copy()
        
        for rr in hh.ascardlist():
            rr.value = interpolate(meta, rr.value)

        prim = pyfits.PrimaryHDU(data=data, header=hh)
        hl = [prim]

        hdulist = pyfits.HDUList(hl)
        return hdulist
