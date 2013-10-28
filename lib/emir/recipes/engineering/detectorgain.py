#
# Copyright 2010-2013 Universidad Complutense de Madrid
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

'''Recipe for the reduction of gain calibration frames.'''

import logging
import math

import numpy
import scipy.stats
import pyfits

import numina.qa
from numina.core import BaseRecipe, Parameter, DataFrame
from numina.core import RecipeError,RecipeRequirements
from numina.core import Product, define_requirements, define_result

from emir.core import RecipeResult
from emir.instrument.detector import CHANNELS, QUADRANTS
from emir.dataproducts import MasterGainMap, MasterRONMap

_logger = logging.getLogger('numina.recipes.emir')

class GainRecipe1Input(RecipeRequirements):
    region = Parameter('channel', 'Region used to compute: (full|quadrant|channel)', 
                       choices=['full','quadrant', 'channel'])
    
class GainRecipe1InputResult(RecipeResult):
    gain = Product(MasterGainMap(None, None, None))
    ron = Product(MasterRONMap(None, None))
    
@define_requirements(GainRecipe1Input)
@define_result(GainRecipe1InputResult)
class GainRecipe1(BaseRecipe):
    '''Detector Gain Recipe.
    
    Recipe to calibrate the detector gain.
    '''

    def __init__(self):
        super(GainRecipe1, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def region(self, reqs):
        mm = reqs['region'].tolower()
        if mm  == 'full':
            return ((slice(0, 2048), slice(0, 2048)))
        elif mm == 'quadrant':
            return QUADRANTS
        elif mm == 'channel':
            return CHANNELS
        else:
            raise ValueError
    
    def run(self, obresult, reqs):

        resets = []
        ramps = []

        for frame in obresult.frames:
            if frame.itype == 'RESET':
                resets.append(frame.label)
                _logger.debug('%s is RESET', frame.label)
            elif frame.itype == 'RAMP':
                ramps.append(frame.label)
                _logger.debug('%s is RAMP', frame.label)
            else:
                raise RecipeError('frame is neither a RAMP nor a RESET')

        channels = self.region(reqs)
        result_gain = numpy.zeros((len(channels), ))
        result_ron = numpy.zeros_like(result_gain)
        
        counts = numpy.zeros((len(ramps), len(channels)))
        variance = numpy.zeros_like(counts)

        last_reset = resets[-1]
        _logger.debug('opening last reset image %s', last_reset)
        last_reset_data = pyfits.getdata(last_reset)

        for i, di in enumerate(ramps):
            with pyfits.open(di, mode='readonly') as fd:
                restdata = fd[0].data - last_reset_data
                for j, channel in enumerate(channels):    
                    c = restdata[channel].mean()
                    _logger.debug('%f counts in channel', c)
                    counts[i, j] = c
                    v = restdata[channel].var(ddof=1)
                    _logger.debug('%f variance in channel', v)
                    variance[i, j] = v

        for j, _ in enumerate(channels):
            slope, intercept, _r_value, _p_value, _std_err = scipy.stats.linregress(counts[:,j], variance[:,j])

            result_gain[j] = 1.0 / slope
            result_ron[j] = math.sqrt(intercept)
        cube = numpy.zeros((2, 2048, 2048))
         
        for gain, var, channel in zip(result_gain, result_ron, channels):
            cube[0][channel] = gain
            cube[1][channel] = var
        
        hdu = pyfits.PrimaryHDU(cube[0])
        hduvar = pyfits.ImageHDU(cube[1])
        hdulist = pyfits.HDUList([hdu, hduvar])

        gain = MasterGainMap(mean=result_gain, var=numpy.array([]), 
                    frame=DataFrame(hdulist))
        ron = MasterRONMap(mean=result_ron, var=numpy.array([]))
        return GainRecipe1InputResult(gain=gain, ron=ron)
        
