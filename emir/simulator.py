#
# Copyright 2008-2011 Universidad Complutense de Madrid
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

from datetime import datetime
from pkgutil import get_data
from StringIO import StringIO
import logging

import pyfits
from pyfits import Header
import numpy
from numpy.random import normal, poisson
from numina.treedict import TreeDict
from numina.instrument.template import interpolate
from numina.recipes import Image

from emir.instrument.detector import EmirDetector

_logger = logging.getLogger('emir.simulator')

class Instrument(object):
    def __init__(self, detector):

        self.detector = detector

        self.meta = TreeDict()
        self.meta['name'] = 'EMIR'
        self.meta['focalstation'] = 'NASMYTH'
        self.meta['filter0'] = 'J'
        self.meta['imagetype'] = ''
        self.meta['detector'] = self.detector.meta
 
    def filter(self, name):
        self.meta['filter0'] = name

    def expose(self, time):
        self.detector.elapse(time)

    def readout(self):
        data = self.detector.readout() 
        return self.meta, data

    def imagetype(self, name):
        self.meta['imagetype'] = name

class Emir(Instrument):
    def __init__(self):
        detector = EmirDetector()
        print detector.meta
        Instrument.__init__(self, detector)

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
        return Header()    
    sfile = StringIO(get_data('emir.instrument', res))
    hh = Header(txtfile=sfile)
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
