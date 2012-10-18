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

'''Tests for RAMP readout mode.'''

from __future__ import division

import unittest

import numpy
import pyfits
from numina.instrument.detector import nIRDetector, Amplifier
from numina.instrument.sources import ThermalBackground

from ..preprocess import axis_ramp, preprocess_ramp
from numina.array._ramp import loopover_ramp_c

class RampReadoutAxisTestCase(unittest.TestCase):
    
    def setUp(self):
        self.emptybp = numpy.zeros((10,), dtype='uint8')
        self.data = numpy.arange(10, dtype='uint16')
        
        self.mask = numpy.zeros((10,), dtype='uint8')
        self.nmap = numpy.zeros((10,), dtype='uint8')
        self.crmask = numpy.zeros((10,), dtype='uint8')

        self.img = numpy.zeros((10,))
        self.var = numpy.zeros((10,))
        self.saturation = 65536
        self.hsize = 5
        self.dt = 1.0
        self.gain = 1.0
        self.ron = 1.0
        self.nsig = 4.0
        self.blank = 0

    def test_saturation0(self):        
        '''Test we count correctly saturated pixels in RAMP mode.'''
        
        MASK_SATURATION = 3 
        MASK_GOOD = 0
    
        # Nno points 
        self.data[:] = 50000 #- 32768
        saturation = 40000 #- 32767
        
        axis_ramp(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.crmask, 
                    saturation, 
                    self.dt, self.gain, self.ron, 
                    self.nsig, 
                    blank=self.blank)

        for nn in self.nmap:
            self.assertEqual(nn, 0)

        for n in self.mask:
            self.assertEqual(n, MASK_SATURATION)
            
        for v in self.var:
            self.assertEqual(v, 0)
            
        for v in self.img:
            self.assertEqual(v, self.blank)
            

    def test_saturation1(self):        
        '''Test we count correctly saturated pixels in RAMP mode.'''
        
        MASK_SATURATION = 3 
        MASK_GOOD = 0
            
        saturation = 50000
        self.data[7:] = saturation 
        
        axis_ramp(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.crmask, 
                    saturation, 
                    self.dt, self.gain, self.ron, 
                    self.nsig, 
                    blank=self.blank)
        
        for nn in self.nmap:
            self.assertEqual(nn, 7)
            
        for n in self.mask:
            self.assertEqual(n, MASK_GOOD)

        for v in self.var:
            self.assertEqual(v, 0.2142857142857143)
            
        for v in self.img:
            self.assertEqual(v, 1)
            
            
        
    def test_badpixel(self):
        '''Test we ignore badpixels in RAMP mode.'''
        self.emptybp[...] = 1

        axis_ramp(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.crmask, 
                    self.saturation, 
                    self.dt, self.gain, self.ron, 
                    self.nsig, 
                    blank=self.blank)
        
        for nn in self.nmap:
            self.assertEqual(nn, 0)
            
        for n in self.mask:
            self.assertEqual(n, 1)
        
        for v in self.var:
            self.assertEqual(v, 0)
            
        for v in self.img:
            self.assertEqual(v, self.blank)
            
            
    def test_results1(self):
        '''Test we obtain correct values in RAMP mode'''
        
        self.data = numpy.arange(10, dtype='int16')

        axis_ramp(self.data, self.emptybp, self.img, self.var, self.nmap, 
                    self.mask, self.crmask, 
                    self.saturation, 
                    self.dt, self.gain, self.ron, 
                    self.nsig, 
                    blank=self.blank)        
                    
        for nn in self.nmap:
            self.assertEqual(nn, 10)
            
        for n in self.mask:
            self.assertEqual(n, 0)
        
        for v in self.var:
            self.assertEqual(v, 0.13454545454545455)
            
        for v in self.img:
            self.assertEqual(v, 1.0)
            
    def test_results2(self):
        '''Test we obtain correct values in RAMP mode'''            
        self.data = numpy.arange(10, dtype='int16')
        
        axis_ramp(self.data, self.emptybp, self.img, self.var, self.nmap, 
                self.mask, self.crmask, 
                self.saturation, 
                self.dt, self.gain, self.ron, 
                self.nsig, 
                blank=self.blank)
                
        for nn in self.nmap:
            self.assertEqual(nn, 10)
            
        for n in self.mask:
            self.assertEqual(n, 0)
        
        for v in self.var:
            self.assertEqual(v, 0.13454545454545455)
            
        for v in self.img:
            self.assertEqual(v, 1.0)
            
    
class RampReadoutFrameTestCase(unittest.TestCase):
    def setUp(self):
        self.frame = None
        
        self.data = numpy.array([[[2343,2454, 2578, 2661,2709, 24311, 24445, 
                                 24405, 24612, 24707]]])
        # self.data is (1,1, 10)
        
        
        hdu = pyfits.PrimaryHDU()
        hdu.data = self.data

        tbcr = 100
        elapsed = 102

        hdu.header.update('readsamp', 10)
        hdu.header.update('exptime', tbcr)
        hdu.header.update('elapsed', elapsed)
        hdu.header.update('readproc', False)
        hdu.header.update('readmode', 'RAMP')

        self.frame = pyfits.HDUList([hdu])
        
    def test_readmode_keyword(self):
        '''Test we raise ValueError if readmode is not RAMP'''
        #preprocess_fowler(frame, saturation=65536, badpixels=None, blank=0):
        
        self.frame[0].header['readmode'] = 'NORAMP'
        
        self.assertRaises(ValueError, preprocess_ramp, self.frame)
        
    def  test_output(self):
        '''Test RAMP output'''
        
        result = preprocess_ramp(self.frame)
        
        self.assertIsInstance(result, pyfits.HDUList)
        
        nex = len(result)
        
        self.assertEqual(nex, 5, 'The number of extensions is not equal to 5')
        
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
        # 5th extension
        ext = result[4]
        self.assertIsInstance(ext, pyfits.ImageHDU)
        self.assertEqual(ext.header.get('extname'), 'CRMASK')
        self.assertEqual(ext.data.dtype, numpy.dtype('uint8'))
        
        
class RampReadoutFrameTestCase2(unittest.TestCase):
    def setUp(self):
        pass
        
    def create_frame(self, source, gain, ron, saturation):
        rows = 6
        columns = 5
        
        chans = [Amplifier((slice(0,rows), slice(0,columns)), gain, ron, saturation)]
        exposure = 10.0
        samples = 10
        
        detector = nIRDetector((rows, columns), chans, dark=0.0, pedestal=0.0)
        detector.readout_time = 0.01
        detector.connect(ThermalBackground(source))
        dt = exposure / samples
        print 'test dt', dt
        print 'ron', ron
        final = numpy.empty((rows, columns, samples), dtype='int32')
        final2 = final[:]
        detector.reset()
        detector.expose(dt)
        final[..., 0] = numpy.random.poisson(dt * source) #detector.readout()
        
        for i in range(1, samples):
            detector.expose(dt)
            final[..., i] = final[..., i -1] + numpy.random.poisson(dt * source) #detector.readout()
            if ron > 0:
                final2[..., i] = numpy.random.normal(final[..., i], ron) / gain
            else:
                final2[..., i] = final[..., i] / gain
            
        final2[final2 < 0] = 0
            
        hdu = pyfits.PrimaryHDU()
        hdu.data = final2

        hdu.header.update('exptime', exposure)
        hdu.header.update('elapsed', detector.time_since_last_reset())
        hdu.header.update('readproc', False)
        hdu.header.update('readmode', 'RAMP')
        hdu.header.update('readsamp', samples)
        hdu.header.update('gain', gain)
        hdu.header.update('ron', ron)

        frame = pyfits.HDUList([hdu])
        return frame 
              
    def test1(self):
        '''Test RAMP output.'''
        
        gain = 1.0
        ron = 0.0
        source = 1000
        saturation = 65000
        
        frame = self.create_frame(source, gain, ron, saturation)
        
        result = preprocess_ramp(frame, saturation=saturation, 
                                 gain=gain, ron=ron)
        
        self.assertIsInstance(result, pyfits.HDUList)
        
        nex = len(result)
        
        print result[0].data
        
        print result[0].data.mean()
        print result[1].data.mean()
        self.assertEqual(nex, 5, 'The number of extensions is not equal to 5')
        
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
        # 5th extension
        ext = result[4]
        self.assertIsInstance(ext, pyfits.ImageHDU)
        self.assertEqual(ext.header.get('extname'), 'CRMASK')
        self.assertEqual(ext.data.dtype, numpy.dtype('uint8'))
        