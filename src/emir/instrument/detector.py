#
# Copyright 2008-2013 Universidad Complutense de Madrid
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

import numpy

from numina.treedict import TreeDict
from numina.instrument.detector import nIRDetector, Channel, DAS
from numina.instrument.detector import SingleReadoutMode
from numina.instrument.detector import CdsReadoutMode
from numina.instrument.detector import RampReadoutMode
from numina.instrument.detector import FowlerReadoutMode

from .channels import CHANNELS, CHANNELS_2, CHANNELS_READOUT
from .channels import CHANNELS_READOUT_2, QUADRANTS

class EMIR_DAS(DAS):
    def __init__(self, detector):
        super(EMIR_DAS, self).__init__(detector)
                
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

class EMIR_Detector_1(nIRDetector):
    def __init__(self, dark=0.25, flat=1.0, bad_pixel_mask=None):
        gain = 3.02
        ron = 6.1 # Electrons
        wdepth = 55292
        pedestal = 2
        resetval = 2
        resetnoise = 0.0
        channels = [Channel((slice(0,2048), slice(0,2048)), gain, ron, wdepth)]
        super(EMIR_Detector_1, self).__init__((2048, 2048), channels, 
                dark=dark, pedestal=pedestal, flat=flat,
                resetval=resetval, resetnoise=resetnoise, 
                bad_pixel_mask=bad_pixel_mask)

class EMIR_Detector_4(nIRDetector):
    def __init__(self, dark=0.25, flat=1.0, bad_pixel_mask=None):
        gain = [3.02, 2.98, 3.00, 2.91]
        ron = [6.1, 5.9, 5.9, 6.2] # Electrons
        wdepth = [55292, 56000, 56000, 56000]
        pedestal = 2
        resetval = 2
        resetnoise = 0.0

        channels = [Channel(*vs) for vs in zip(QUADRANTS, gain, ron, wdepth)]
        super(EMIR_Detector_4, self).__init__((2048, 2048), channels, 
                dark=dark, pedestal=pedestal, flat=flat,
                resetval=resetval, resetnoise=resetnoise, 
                bad_pixel_mask=bad_pixel_mask)

class EMIR_Detector_32(nIRDetector):
    def __init__(self, dark=0.25, flat=1.0, bad_pixel_mask=None):
        # ADU, per AMP
        ron = [3.14, 3.060, 3.09, 3.08, 3.06, 3.20, 3.13, 3.10, 2.98, 2.98, 
               2.95, 3.06, 3.050, 3.01, 3.05, 3.20, 2.96, 3.0, 3.0, 2.99, 3.14, 
               2.99, 3.06, 3.05, 3.08, 3.06, 3.01, 3.02, 3.07, 2.99, 3.03, 3.06]

        # e-/ADU per AMP
        gain = [2.83, 2.84, 2.82, 2.92, 2.98, 2.71, 3.03, 2.92, 3.04, 2.79, 2.72, 
                3.07, 3.03, 2.91, 3.16, 3.22, 2.93, 2.88, 2.99, 3.24, 3.25, 3.12, 
                3.29, 3.03, 3.04, 3.35, 3.11, 3.25, 3.29, 3.17, 2.93, 3.05]

        # ADU per AMP
        wdepth = [42353.1, 42148.3, 42125.5, 42057.9, 41914.1, 42080.2, 42350.3, 
                  41830.3, 41905.3, 42027.9, 41589.5, 41712.7, 41404.9, 41068.5, 
                  40384.9, 40128.1, 41401.4, 41696.5, 41461.1, 41233.2, 41351.0, 
                  41803.7, 41450.2, 41306.2, 41609.4, 41414.1, 41324.5, 41691.1, 
                  41360.0, 41551.2, 41618.6, 41553.5]
        
        pedestal = 2
        resetval = 2
        resetnoise = 0.0

        channels = [Channel(*vs) for vs in zip(CHANNELS, gain, ron, wdepth)]
        super(EMIR_Detector_32, self).__init__((2048, 2048), channels, 
                dark=dark, pedestal=pedestal, flat=flat,
                resetval=resetval, resetnoise=resetnoise, 
                bad_pixel_mask=bad_pixel_mask)

EMIR_Detector = EMIR_Detector_32
