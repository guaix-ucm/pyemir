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

from numina.instrument.detector import nIRDetector, Channel, DAS
from numina.instrument.detector import SingleReadoutMode
from numina.instrument.detector import CdsReadoutMode
from numina.instrument.detector import RampReadoutMode
from numina.instrument.detector import FowlerReadoutMode

from .channels import CHANNELS_3
from .channels import CHANNELS
from .channels import QUADRANTS

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
        ron = 2.1 # ADU
        pedestal = 5362
        wdepth = 55292
        saturation = 57000
        channels = [Channel((slice(0,2048), slice(0,2048)), gain, ron, pedestal, wdepth, saturation)]
        super(EMIR_Detector_1, self).__init__((2048, 2048), channels, 
                dark=dark, flat=flat,
                bad_pixel_mask=bad_pixel_mask)

class EMIR_Detector_4(nIRDetector):
    def __init__(self, dark=0.25, flat=1.0, bad_pixel_mask=None):
        gain = [3.02, 2.98, 3.00, 2.91]
        ron = [2.1, 1.9, 1.9, 2.2] # Electrons
        pedestal = [5362, 5362, 5362, 5362]
        wdepth = [55292, 56000, 56000, 56000]
        saturation = [57292, 57000, 57000, 57000]

        channels = [Channel(*vs) for vs in zip(QUADRANTS, gain, ron, pedestal, 
        wdepth, saturation)]
        super(EMIR_Detector_4, self).__init__((2048, 2048), channels, 
                dark=dark, flat=flat,
                bad_pixel_mask=bad_pixel_mask)

class EMIR_Detector_32(nIRDetector):
    def __init__(self, dark=0.25, flat=1.0, bad_pixel_mask=None):
        # from Carlos Gonzalez Thesis, page 147
        # ADU
        ron = [2.03, 1.98, 1.96, 1.95, 1.99, 1.95, 1.97, 1.95,
               1.94, 1.96, 1.94, 1.92, 2.03, 1.93, 1.96, 1.98,
               2.01, 1.99, 1.97, 1.98, 1.99, 1.97, 1.97, 2.00,
               2.02, 2.02, 2.05, 1.98, 2.00, 2.02, 2.01, 2.03
              ]

        # from Carlos Gonzalez Thesis, page 127
        # e-/ADU
        gain = [3.08, 2.90, 2.68, 3.12, 2.63, 3.10, 3.00, 3.02,
                3.18, 3.11, 3.09, 3.19, 3.11, 2.99, 2.60, 3.02,
                2.99, 3.16, 3.18, 3.11, 3.17, 3.07, 3.12, 3.02,
                2.92, 3.07, 2.90, 2.91, 2.95, 3.00, 3.01, 3.01
                ]
        
        # from Carlos Gonzalez Thesis, page 130
        # ADU
        pedestal = [3776, 4428, 4071, 4161, 5540, 5623, 5819, 5933,
                    6245, 6151, 5746, 5477, 4976, 4937, 4966, 4576,
                    5605, 5525, 5532, 5643, 5505, 5371, 5535, 4985,
                    5651, 5707, 5604, 5260, 4696, 5053, 5032, 5155
                    ]

        # from Carlos Gonzalez Thesis, page 132
        # ADU
        wdepth = [54110, 54796, 55205, 55349, 55505, 54990, 55360, 55276,
                  55071, 55585, 55419, 55126, 55284, 55884, 55523, 55416,
                  55587, 55356, 55273, 55705, 55203, 55356, 55275, 55235,
                  55123, 55549, 54845, 54895, 55207, 55340, 55556, 56009
        ]

        # from Carlos Gonzalez Thesis, page 132
        # ADU
        saturation = [56811, 57011, 57042, 57057, 57178, 56738, 57171, 57094,
                      56866, 57215, 57022, 56761, 56983, 57655, 57329, 57212,
                      57299, 57047, 56943, 57377, 56878, 57035, 57048, 57072,
                      57224, 57829, 57407, 57178, 57398, 57496, 57564, 57938
                      ]

        channels = [Channel(*vs) for vs in zip(CHANNELS_3, gain, ron, pedestal, wdepth, saturation)]
        super(EMIR_Detector_32, self).__init__((2048, 2048), channels, 
                dark=dark,flat=flat,
                bad_pixel_mask=bad_pixel_mask)

EMIR_Detector = EMIR_Detector_32
