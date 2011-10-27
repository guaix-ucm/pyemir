#
# Copyright 2010-2011 Sergio Pascual
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
from numina.recipes import RecipeBase
#from emir.dataproducts import create_result
from emir.recipes import EmirRecipeMixin
from numina.recipes.registry import ProxyPath
from numina.recipes.registry import Schema
from numina.image import DiskImage
from emir.instrument.detector import CHANNELS, QUADRANTS
from emir.dataproducts import create_result

_logger = logging.getLogger("emir.recipes")


def multimap(fun, ll):
    if hasattr(ll, "__iter__"):
        return [multimap(fun, i) for i in ll]
    return fun(ll)

class Recipe1(RecipeBase, EmirRecipeMixin):
    '''Detector Gain Recipe.
    
    Recipe to calibrate the detector gain.
    '''
    required_parameters = [
        Schema('resets', ProxyPath('/observing_block/result/resets'), 
               'A list of paths to reset images'),
        Schema('ramps', ProxyPath('/observing_block/result/ramps'), 
               'A list of ramps'),
        Schema('region', 'channel', 'Region used to compute: (full|quadrant|channel)')
    ]
    capabilities = ['detector_gain']
    
    
    def __init__(self, param, runinfo):
        super(Recipe1, self).__init__(param, runinfo)
        _logger.info('building')

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

    def run(self):
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
        
class Recipe2(RecipeBase, EmirRecipeMixin):
    '''Detector Gain Recipe.
    
    Recipe to calibrate the detector gain.
    '''
    required_parameters = [
        Schema('value', 0, 'Number')
    ]
    capabilities = ['detector_gain']
    
    
    def __init__(self, param, runinfo):
        super(Recipe2, self).__init__(param, runinfo)

    def run(self):  
        
        return {'qa': numina.qa.UNKNOWN, 'gain': {'mean': self.parameters['value'], 
                                                  'var': 0},
                                          }
        
        
class Recipe3(RecipeBase, EmirRecipeMixin):
    '''Detector Gain Recipe.
    
    Recipe to calibrate the detector gain.
    '''
    required_parameters = [
        Schema('value', 0, 'Number')
    ]
    capabilities = ['detector_gain']
    
    
    def __init__(self, param, runinfo):
        super(Recipe3, self).__init__(param, runinfo)

    def run(self):  
        
        return {'qa': numina.qa.UNKNOWN, 'gain': {'mean': self.parameters['value'], 
                                                  'var': 0},
                                          }
        
    
if __name__ == '__main__':
    import os
    
    from numina.compatibility import json
    from numina.user import main
    from numina.jsonserializer import to_json

    logging.basicConfig(level=logging.DEBUG)
    _logger.setLevel(logging.DEBUG)
        
    pv = {'recipes': {'default': {'parameters': {
                                    'region': 'channel'
                                },
                     'run': {'repeat': 1,
                             },   
                        },
                      'emir.recipes.detector_gain.Recipe1': {
                        'parameters': {
                           'region': 'full'
                    },                                                            
                                                            },
                      'emir.recipes.detector_gain.Recipe2': {
                        'parameters': {
                           'value': 0
                    },
                        'run': {'repeat': 2,
                             },   
                        },                 
                      },
                      
          'observing_block': {'instrument': 'emir',
                       'mode': 'detector_gain',
                       'id': 1,
                       'result': {
                                  'resets': ['apr21_0046.fits', 
                                             'apr21_0047.fits', 
                                             'apr21_0048.fits',
                                             'apr21_0049.fits', 
                                             'apr21_0050.fits'
                                             ],
                                   'ramps': [
                                             ['r0000.fits', 'r0001.fits','r0002.fits', 'r0003.fits',
                                               'r0004.fits', 'r0005.fits', 'r0006.fits', 'r0007.fits'],
                                            ['r0008.fits', 'r0009.fits','r0010.fits', 'r0011.fits',
                                               'r0012.fits', 'r0013.fits', 'r0014.fits', 'r0015.fits'],
                                               
                                             ]
                            },
                       },                     
    }
        
    os.chdir('/home/spr/Datos/emir/test5')
    
    ff = open('config.txt', 'w+')
    try:
        json.dump(pv, ff, default=to_json, encoding='utf-8', indent=2)
    finally:
        ff.close()

    main(['-d', '--resultsdir', '/home/spr', '--run', 'config.txt'])
