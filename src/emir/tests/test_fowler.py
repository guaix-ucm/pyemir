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

'''Tests for Fowler readout mode.'''

import unittest

import numpy
import pyfits

from ..preprocess import axis_fowler, preprocess_fowler

class FowlerReadoutAxisTestCase(unittest.TestCase):
    
    def setUp(self):
        self.emptybp = numpy.zeros((10,), dtype='uint8')
        self.data = numpy.zeros((10,), dtype='uint16')
        
        self.mask = numpy.zeros((10,), dtype='uint8')
        self.nmap = numpy.zeros((10,), dtype='uint8')

        self.img = numpy.zeros((10,))
        self.var = numpy.zeros((10,))
        self.saturation = 65536
        self.hsize = 5
        self.blank = 0
    
    def test_saturation(self):        
        '''Test we count correctly saturated pixels in Fowler mode.'''
        
        MASK_SATURATION = 3 
        MASK_GOOD = 0
    
        # 5 to 9 are saturated, no points 
        self.data[5:] = 50000 #- 32768
        saturation = 40000 #- 32767
        
        axis_fowler(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.hsize, 
                    saturation, blank=self.blank)

        for nn in self.nmap:
            self.assertEqual(nn, 0)

        for n in self.mask:
            self.assertEqual(n, MASK_SATURATION)
            
        for v in self.var:
            self.assertEqual(v, 0)
            
        for v in self.img:
            self.assertEqual(v, self.blank)
            
            
        # 5 and 6 are OK, 7,8 and 9 are saturated, 2 points
        val = 400
        
        self.data[7:] = 50000
        self.data[:5] = 0
        self.data[5:7] = 400 
        
        axis_fowler(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.hsize, 
                    saturation, blank=self.blank)            
            
            
        for nn in self.nmap:
            self.assertEqual(nn, 2)
            
        for n in self.mask:
            self.assertEqual(n, MASK_GOOD)

        for v in self.var:
            self.assertEqual(v, 0)
            
        for v in self.img:
            self.assertEqual(v, val)
            
        
    def test_badpixel(self):
        '''Test we ignore badpixels in Fowler mode.'''
        self.emptybp[...] = 1

        axis_fowler(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.hsize, 
                    self.saturation, blank=self.blank)
        
        for nn in self.nmap:
            self.assertEqual(nn, 0)
            
        for n in self.mask:
            self.assertEqual(n, 1)
        
        for v in self.var:
            self.assertEqual(v, 0)
            
        for v in self.img:
            self.assertEqual(v, self.blank)
            
            
    def test_results(self):
        '''Test we obtain correct values in Fowler mode'''
        
        self.data = numpy.array([3,4,5,6,7, 403, 404, 405, 406, 407])
        
        axis_fowler(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.hsize, 
                    self.saturation, blank=self.blank)
            
        
        for nn in self.nmap:
            self.assertEqual(nn, 5)
            
        for n in self.mask:
            self.assertEqual(n, 0)
        
        for v in self.var:
            self.assertEqual(v, 0)
            
        for v in self.img:
            self.assertEqual(v, 400)
            
            
        self.data = numpy.array([2343,2454, 2578, 2661,2709, 24311, 24445, 
                                 24405, 24612, 24707])
        
        axis_fowler(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.hsize, 
                    self.saturation, blank=self.blank)
            
        
        for nn in self.nmap:
            self.assertEqual(nn, 5)
            
        for n in self.mask:
            self.assertEqual(n, 0)
        
        for v in self.var:
            self.assertEqual(v, 775.76)
            
        for v in self.img:
            self.assertEqual(v, 21947.0)
            
class FowlerReadoutFrameTestCase(unittest.TestCase):
    def setUp(self):
        self.frame = None
        
        self.data = numpy.array([[[2343,2454, 2578, 2661,2709, 24311, 24445, 
                                 24405, 24612, 24707]]])
        # self.data is (1,1, 10)
        
        
        hdu = pyfits.PrimaryHDU()
        hdu.data = self.data

        tbcr = 100
        elapsed = 101

        hdu.header.update('exptime', tbcr)
        hdu.header.update('elapsed', elapsed)
        hdu.header.update('readproc', False)
        hdu.header.update('readmode', 'FOWLER')

        self.frame = pyfits.HDUList([hdu])
        
    def test_readmode_keyword(self):
        '''Test we raise ValueError if readmode is not FOWLER'''
        #preprocess_fowler(frame, saturation=65536, badpixels=None, blank=0):
        
        self.frame[0].header['readmode'] = 'NOFOWLER'
        
        self.assertRaises(ValueError, preprocess_fowler, self.frame)
        
    def  test_output(self):
        '''Test Fowler output'''
        
        result = preprocess_fowler(self.frame)
        
        self.assertIsInstance(result, pyfits.HDUList)
        
        nex = len(result)
        
        self.assertEqual(nex, 4, 'The number of extensions is not equal to 4')
        
        # first extension
        ext = result[0]
        self.assertTrue(ext.header['readproc'])
        
        # second extension
        ext = result[1]
        self.assertIsInstance(ext, pyfits.ImageHDU)
        self.assertEqual(ext.header.get('extname'), 'VARIANCE')
        # 3rd extension
        ext = result[2]
        self.assertIsInstance(ext, pyfits.ImageHDU)
        self.assertEqual(ext.header.get('extname'), 'MAP')
        self.assertEqual(ext.data.dtype, numpy.dtype('uint8'))
        # 4th extension
        ext = result[3]
        self.assertIsInstance(ext, pyfits.ImageHDU)
        self.assertEqual(ext.header.get('extname'), 'MASK')
        self.assertEqual(ext.data.dtype, numpy.dtype('uint8'))
        
        
        
        
        
        