#
# Copyright 2010-2012 Universidad Complutense de Madrid
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
import os.path

import numpy
import scipy.stats

import numina.qa
from numina.recipes import RecipeBase, provides, Parameter
#from emir.dataproducts import create_result

from emir.instrument.detector import CHANNELS, QUADRANTS

from ...dataproducts import create_result

from ...dataproducts import MasterGainMap, MasterRONMap

_logger = logging.getLogger("emir.recipes")


def multimap(fun, ll):
    if hasattr(ll, "__iter__"):
        return [multimap(fun, i) for i in ll]
    return fun(ll)

@provides(MasterGainMap, MasterRONMap)
class Recipe1(RecipeBase):
    '''Detector Gain Recipe.
    
    Recipe to calibrate the detector gain.
    '''
    __requires__ = [
#        Schema('resets', ProxyPath('/observing_block/result/resets'), 
#               'A list of paths to reset images'),
#        Schema('ramps', ProxyPath('/observing_block/result/ramps'), 
#               'A list of ramps'),
        Parameter('region', 'channel', 'Region used to compute: (full|quadrant|channel)')
    ]

    def __init__(self):
        super(Recipe1, self).__init__(
            author="Sergio Pascual <sergiopr@fis.ucm.es>",
            version="0.1.0"
        )

    def run(self, obresult):
        return {'products': [MasterGainMap(), MasterRONMap()]}
  
    
    def region_full(self):
        return [(slice(0, 2048), slice(0, 2048))]
    
    def region_quadrant(self):
        return QUADRANTS
    
    def region_channel(self):
        return CHANNELS

    def region(self):
        fun = getattr(self, 'region_%s' %  self.parameters['region'])  
        return fun()
    
    def setup(self):
        # Sanity check, check: all images belong to the same detector mode
        
        self.parameters['resets'] = map(lambda x: DiskImage(os.path.abspath(x)), 
                                             self.parameters['resets'])
        self.parameters['ramps'] = multimap(lambda x: DiskImage(os.path.abspath(x)), 
                                             self.parameters['ramps'])

    def old_run(self):
        channels = self.region()
        ramps = self.parameters['ramps']
        result_gain = numpy.zeros((len(ramps), len(channels)))
        result_ron = numpy.zeros_like(result_gain)
        
        for ir, ramp in enumerate(ramps):
            counts = numpy.zeros((len(ramp), len(channels)))
            variance = numpy.zeros_like(counts)
            for i, di in enumerate(ramp):
                di.open(mode='readonly')
                try:
                    for j, channel in enumerate(channels):    
                        counts[i][j] = di.data[channel].mean()
                        variance[i][j] = di.data[channel].var(ddof=1)
                finally:
                    di.close()

            for j, _ in enumerate(channels):
                ig, ron,_,_,_ = scipy.stats.linregress(counts[:,j], variance[:,j])

                result_gain[ir][j] = 1.0 / ig
                result_ron[ir][j] = ron
        _logger.info('running')
        gch_mean = result_gain.mean(axis=0)
        gch_var = result_gain.var(axis=0, ddof=1)
        rch_mean = result_ron.mean(axis=0)
        rch_var = result_ron.var(axis=0, ddof=1)
        
        cube = numpy.zeros((2, 2048, 2048))
         
        for gain, var, channel in zip(gch_mean, gch_var, channels):
            cube[0][channel] = gain
            cube[1][channel] = var
        
        result = create_result(cube[0], variance=cube[1])                                    
        
        return {'qa': numina.qa.UNKNOWN, 'gain': {'mean': list(gch_mean.flat), 
                                                  'var': list(gch_var.flat),
                                                  'image': result
                                                  },
                                          'ron': {'mean': list(rch_mean.flat), 
                                                  'var': list(rch_var.flat)},
                                                  }
