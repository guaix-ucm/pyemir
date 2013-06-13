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
from numina.instrument.detector import nIRDetector, Channel
from numina.instrument.sources import ThermalBackground

from ..preprocess import preprocess_ramp

class RampReadoutFrameTestCase(unittest.TestCase):
    def setUp(self):
        self.frame = None
        
        self.data = numpy.array([[[2343,2454, 2578, 2661,2709, 24311, 24445, 
                                 24405, 24612, 24707]]])
        # self.data is (1,1, 10)
        self.gain = 1.0
        self.ron = 1.0
        
        
        hdu = pyfits.PrimaryHDU()
        hdu.data = self.data

        tbcr = 100
        elapsed = 102

        hdu.header.update('readsamp', 10)
        hdu.header.update('exptime', tbcr)
        hdu.header.update('exposed', tbcr)
        hdu.header.update('elapsed', elapsed)
        hdu.header.update('readproc', False)
        hdu.header.update('readmode', 'RAMP')

        self.frame = pyfits.HDUList([hdu])
        
    def test_readmode_keyword(self):
        '''Test we raise ValueError if readmode is not RAMP'''
        #preprocess_fowler(frame, saturation=65536, badpixels=None, blank=0):
        
        self.frame[0].header['readmode'] = 'NORAMP'
        
        self.assertRaises(ValueError, preprocess_ramp, self.frame, self.gain, self.ron)
        
    def  test_output(self):
        '''Test RAMP output'''
        
        result = preprocess_ramp(self.frame, self.gain, self.ron)
        
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
        
    def create_frame(self, source, gain, ron, pedestal, wdepth, saturation):
        rows = 6
        columns = 5
        
        chans = [Channel((slice(0,rows), slice(0,columns)), gain, ron, pedestal, wdepth, saturation)]
        exposure = 10.0
        samples = 10
        
        detector = nIRDetector((rows, columns), chans, dark=0.0)
        detector.readout_time = 0.01
        detector.connect(ThermalBackground(source))
        dt = exposure / samples
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
        hdu.header.update('exposed', exposure)
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
        pedestal = 65000
        wdepth = 65000
        saturation = 65000
        
        frame = self.create_frame(source, gain, ron, pedestal, wdepth, saturation)
        
        result = preprocess_ramp(frame, gain, ron, saturation=saturation)
        
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


