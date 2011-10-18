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

import itertools as ito

import numpy # pylint: disable-msgs=E1101

from numina import braid
from numina.instrument.detector import Detector

# Classes are new style
__metaclass__ = type

def _channel_gen1(beg, end, step):
    return ito.imap(lambda x: (x, x + step), xrange(beg, end, step))

def _channel_gen2(beg, end, step):
    return ito.imap(lambda x: (x - step, x), xrange(beg, end, -step))

def _ch1():
    return ito.izip(ito.repeat(slice(1024, 2048)), ito.starmap(slice, _channel_gen2(1024, 0, 128)))

def _ch2():  
    return ito.izip(ito.starmap(slice, _channel_gen2(1024, 0, 128)), ito.repeat(slice(0, 1024)))

def _ch3():
    return ito.izip(ito.repeat(slice(0, 1024)), ito.starmap(slice, _channel_gen1(1024, 2048, 128)))

def _ch4():
    return ito.izip(ito.starmap(slice, _channel_gen1(1024, 2048, 128)), ito.repeat(slice(1024, 2048)))

# Channels are listed per quadrant and then in fast readout order
CHANNELS = list(ito.chain(_ch1(), _ch2(), _ch3(), _ch4()))

# Channels in read out order
CHANNELS_READOUT = list(braid(_ch1(), _ch2(), _ch3(), _ch4()))

# Quadrants are listed starting at left-top and counter-clockwise then
QUADRANTS = [(slice(1024, 2048), slice(0, 1024)),
             (slice(0, 1024), slice(0, 1024)),
             (slice(0, 1024), slice(1024, 2048)),
             (slice(1024, 2048), slice(1024, 2048))
             ]

class ReadoutMode:
    def __init__(self, mode, scheme, reads, repeat):
        self.mode = mode
        self.scheme = scheme
        self.reads = reads
        self.repeat = repeat
        
class SingleReadoutMode(ReadoutMode):
    def __init__(self, repeat=1):
        ReadoutMode.__init__(self, 'single', 'perline', 1, repeat)
        
    def events(self, exposure):
        return [exposure]
    
    def process(self, images, events):
        return images[0]
    
class CdsReadoutMode(ReadoutMode):
    '''Correlated double sampling readout mode.'''
    def __init__(self, repeat=1):
        ReadoutMode.__init__(self, 'CDS', 'perline', 1, repeat)
        
    def events(self, exposure):
        '''
        >>> cds_rm = CdsReadoutMode()
        >>> cds_rm.events(10.0)
        [0.0, 10.0]
        '''
        return [0.0, exposure]
    
    def process(self, images, events):
        return images[1] - images[0]
    
class FowlerReadoutMode(ReadoutMode):
    '''Fowler sampling readout mode.'''
    def __init__(self, reads, repeat=1, readout_time=0.0):
        ReadoutMode.__init__(self, 'Fowler', 'perline', reads, repeat)
        self.readout_time = readout_time
        
    def events(self, exposure):
        '''
        
        >>> frm = FowlerReadoutMode(reads=3, readout_time=0.8)
        >>> frm.events(10.0)
        [0.0, 0.8, 1.6, 10.0, 10.8, 11.6]
        '''
        
        dt = self.readout_time
        vreads = [i * dt for i in range(self.reads)]
        vreads += [t + exposure for t in vreads]
        return vreads
    
    def process(self, images, events):
        # Subtracting correlated reads
        nsamples = len(images) / 2
        # nsamples has to be odd
        reduced = numpy.array([images[nsamples + i] - images[i] 
                               for i in range(nsamples)])
        # Final mean
        return reduced.mean(axis=0)

class RampReadoutMode(ReadoutMode):
    '''"Up the ramp" sampling readout mode.'''
    def __init__(self, reads, repeat=1):
        ReadoutMode.__init__(self, 'Ramp', 'perline', reads, repeat)
        
    def events(self, exposure):
        dt = exposure / (self.reads - 1.)
        return [dt * i for i in range(self.reads)]
    
    def process(self, images, events):
        
        def slope(y, xcenter, varx, time):
            return ((y - y.mean()) * xcenter).sum() / varx * time
        
        events = numpy.asarray(events)
        images = numpy.asarray(images)
        meanx = events.mean()
        sxx = events.var() * events.shape[0]
        xcenter = events - meanx
        images = numpy.dstack(images)
        return numpy.apply_along_axis(slope, 2, images, xcenter, sxx, events[ - 1])

class Hawaii2Detector(Detector):
    '''Hawaii2 detector.'''
    
    AMP1 = QUADRANTS # 1 amplifier per quadrant
    AMP8 = CHANNELS # 8 amplifiers per quadrant
    
    shape = (2048, 2048)
    amplifiers = CHANNELS
    
    def __init__(self, gain=1.0, ron=0.0, dark=1.0, well=65535,
                 pedestal=200., flat=1.0, resetval=0, resetnoise=0.0,
                 mode='8'):
        '''
            :parameter gain: gain in e-/ADU
            :parameter ron: ron in ADU
            :parameter dark: dark current in e-/s
            :parameter well: well depth in ADUs 
        '''
        super(Hawaii2Detector, self).__init__(self.shape, gain, ron, dark, well,
                                           pedestal, flat, resetval, resetnoise)
        
        if mode not in ['1', '8']:
            raise ValueError('mode must be "1" or "8"')
        
        self.mode = mode
        # Amplifier region
        self.amplifiers = self.AMP1 if mode == '1' else self.AMP8
        
        # Gain and RON per amplifier
        self._ron = numpy.asarray(ron)
        self._gain = numpy.asarray(gain)
        
        self.ronmode = SingleReadoutMode()
        
        self.events = None
        
        self._exposure = 0
        
    def read(self, time=None, source=None):
        '''Read the detector.'''
        if time is not None:
            self.elapse(time, source)
        self._time += self.readout_time
        result = self._detector.copy()
        result[result < 0] = 0
        
        # Gain and RON per amplifier        
        ampgain = ito.cycle(self._gain.flat)
        ampron = ito.cycle(self._ron.flat)
        
        for amp, gain, ron in zip(self.amplifiers, ampgain, ampron):
            data = result[amp]
            data /= gain            
            # Readout noise
            data += numpy.random.standard_normal(data.shape) * ron
        
        result += self._pedestal
        # result[result > self._well] = self._well
        return result.astype(self.type)        
        
    def configure(self, ronmode):
        self.ronmode = ronmode
    
    def exposure(self, exposure):
        self._exposure = exposure
        self.events = self.ronmode.events(exposure)
    
    def lpath(self, source=None):
        self.reset()
        images = [self.read(t, source) for t in self.events]
        # Process the images according to the mode
        final = self.ronmode.process(images, self.events)
        return final.astype(self.outtype)

    def metadata(self):
        '''Return metadata exported by the EmirDetector.'''
        mtdt = {'EXPOSED':self._exposure, 
                'EXPTIME':self._exposure,
                'ELAPSED':self.time_since_last_reset(),
                'DARKTIME':self.time_since_last_reset(),
                'READMODE':self.ronmode.mode.upper(),
                'READSCHM':self.ronmode.scheme.upper(),
                'READNUM':self.ronmode.reads,
                'READREPT':self.ronmode.repeat}
        return mtdt

class EmirDetector(Hawaii2Detector):
    def __init__(self, flat=1.0):        
        # ADU, per AMP
        ron = [3.14, 3.060, 3.09, 3.08, 3.06, 3.20, 3.13, 3.10, 2.98, 2.98, 
               2.95, 3.06, 3.050, 3.01, 3.05, 3.20, 2.96, 3.0, 3.0, 2.99, 3.14, 
               2.99, 3.06, 3.05, 3.08, 3.06, 3.01, 3.02, 3.07, 2.99, 3.03, 3.06]

        # e-/ADU per AMP
        gain = [2.83, 2.84, 2.82, 2.92, 2.98, 2.71, 3.03, 2.92, 3.04, 2.79, 2.72, 
                3.07, 3.03, 2.91, 3.16, 3.22, 2.93, 2.88, 2.99, 3.24, 3.25, 3.12, 
                3.29, 3.03, 3.04, 3.35, 3.11, 3.25, 3.29, 3.17, 2.93, 3.05]

        # e-/s global median
        dark = 0.298897832632

        # ADU per AMP
        wdepth = [42353.1, 42148.3, 42125.5, 42057.9, 41914.1, 42080.2, 42350.3, 
                  41830.3, 41905.3, 42027.9, 41589.5, 41712.7, 41404.9, 41068.5, 
                  40384.9, 40128.1, 41401.4, 41696.5, 41461.1, 41233.2, 41351.0, 
                  41803.7, 41450.2, 41306.2, 41609.4, 41414.1, 41324.5, 41691.1, 
                  41360.0, 41551.2, 41618.6, 41553.5]

        super(EmirDetector, self).__init__(gain=gain, ron=ron, dark=dark, 
                                           well=wdepth, flat=flat)

class Hawaii1Detector(Detector):
    '''Hawaii1 detector.'''
        
    shape = (1024, 1024)
    amplifiers =  [(slice(512, 1024), slice(0, 512)),
             (slice(0, 512), slice(0, 512)),
             (slice(0, 512), slice(512, 1024)),
             (slice(512, 1024), slice(512, 1024))
             ]

