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

import itertools as ito

import numpy # pylint: disable-msgs=E1101

from numina.extraiter import braid
from numina.treedict import TreeDict
from numina.instrument.detector import nIRDetector, Amplifier, Das
from numina.instrument.detector import SingleReadoutMode
from numina.instrument.detector import CdsReadoutMode
from numina.instrument.detector import RampReadoutMode
from numina.instrument.detector import FowlerReadoutMode

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

class EmirDas(Das):
    def __init__(self, detector):
        Das.__init__(self, detector)
                
    def readout_mode_single(self, repeat=1):
        self.readmode(SingleReadoutMode(repeat=repeat)) 
    
    def readout_mode_cds(self, repeat=1):
        mode = CdsReadoutMode(repeat=repeat)
        self.readmode(mode)
    
    def readout_mode_fowler(self, reads, repeat=1, readout_time=0.0):
        mode = FowlerReadoutMode(reads, repeat=repeat, 
                                 readout_time=readout_time)
        self.readmode(mode)    

    def readout_mode_ramp(self, reads, repeat=1):
        mode = RampReadoutMode(reads, repeat=repeat)
        self.readmode(mode)

class Hawaii2Detector(nIRDetector):
    '''Hawaii2 detector.'''
    
    AMP1 = QUADRANTS # 1 amplifier per quadrant
    AMP8 = CHANNELS # 8 amplifiers per quadrant
    
    shape = (2048, 2048)
    amplifiers = CHANNELS
    
    def __init__(self, gain=1.0, ron=0.0, dark=1.0, well=65535,
                 pedestal=200., flat=1.0, resetval=0, resetnoise=0.0,
                 ampmode='8'):
        '''
            :parameter gain: gain in e-/ADU
            :parameter ron: ron in ADU
            :parameter dark: dark current in e-/s
            :parameter well: well depth in ADUs 
        '''
        
        
        if ampmode not in ['1', '8']:
            raise ValueError('ampmode must be "1" or "8"')
        
        self.ampmode = ampmode
        # Amplifier region
        self.ampgeom = self.AMP1 if ampmode == '1' else self.AMP8
        
        ampgain = ito.cycle(numpy.asarray(gain).flat)
        ampron = ito.cycle(numpy.asarray(ron).flat)
        ampwell = ito.cycle(numpy.asarray(ron).flat)
        amps = [Amplifier(geom, gain, ron, well) for geom, gain, ron, well in zip(self.amplifiers, ampgain, ampron, ampwell)]
        
        nIRDetector.__init__(self, self.shape, amps, dark, pedestal, flat, resetval, resetnoise)
        
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
        
        self.meta = TreeDict()
        self.meta['gain'] = 2.8
        self.meta['readnoise'] = 3.0

        Hawaii2Detector.__init__(self, gain=gain, ron=ron, dark=dark, 
                                           well=wdepth, flat=flat)

if __name__ == '__main__':
    
    
    det = EmirDetector()
