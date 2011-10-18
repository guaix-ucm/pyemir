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

'''Recipe for the reordering of frames.'''

import logging
import os.path
import itertools as ito

import numpy # pylint: disable-msgs=E1101


import numina.qa
from numina.recipes import RecipeBase
from emir.recipes import EmirRecipeMixin
from numina.recipes.registry import ProxyPath
from numina.recipes.registry import Schema
from numina.image import DiskImage
from emir.instrument.detector import CHANNELS_READOUT
from emir.dataproducts import create_result

_logger = logging.getLogger("emir.recipes")

class Recipe(RecipeBase, EmirRecipeMixin):
    '''Reordering Recipe.
    
    Recipe to reorder images created by the detector.
    '''
    required_parameters = [
        Schema('images', ProxyPath('/observing_block/result/images'), 
               'A list of images created by the detector'),
    ]
    capabilities = ['detector_reorder']
    
    
    def __init__(self, param, runinfo):
        super(Recipe, self).__init__(param, runinfo)
        
    def setup(self):
        self.parameters['images'] = [DiskImage(os.path.abspath(path)) 
                  for path in self.parameters['images']]

    def run(self):
        
        images = self.parameters['images']
        
        results = []

        for img in images:
            # Using numpy memmap instead of pyfits
            # Images are a in a format unrecognized by pyfits
            _logger.debug('processing %s', img.filename)
            f = numpy.memmap(img.filename, 
                                  dtype='uint16', mode='r', offset=36 * 80)            
            try:                                
                f.shape = (1024, 4096)
                rr = numpy.zeros((2048, 2048), dtype='int16')
                for idx, channel, conv in ito.izip(ito.count(0), 
                                                   CHANNELS_READOUT,
                                                   ito.cycle([lambda x:x, lambda x:x.T])
                                                   ):
                    rr[channel] = conv(f[:,slice(128 * idx, 128 * (idx + 1))])

                basename = os.path.basename(img.filename)
                primary_headers = {'FILENAME': basename}
                result = create_result(rr, 
                                       headers=primary_headers)
                result.writeto(basename)
                
                result.close()                
                results.append(DiskImage(os.path.abspath(basename)))
            finally:
                del f
                img.close()
        
        return {'qa': numina.qa.UNKNOWN, 'images': results}
                                                  
    
if __name__ == '__main__':
    import uuid
    import glob
    
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
                      },                      
          'observing_block': {'instrument': 'emir',
                       'mode': 'detector_reorder',
                       'id': 1,
                       'result': {
                                  'images': [], 

                            },
                       },                     
    }
    os.chdir('/home/spr/Datos/emir/test6/data')
    
    pv['observing_block']['result']['images'] = glob.glob('*serie*.fits')
    #pv['observing_block']['result']['images'] = ['pozo_BG33_VR05_OFF11_serie.0054.fits']
    os.chdir('/home/spr/Datos/emir/test6')
    
    # Creating base directory for storing results
    uuidstr = str(uuid.uuid1()) 
    basedir = os.path.join(os.getcwd(), uuidstr)
    os.mkdir(basedir)
    
    ff = open('config.txt', 'w+')
    try:
        json.dump(pv, ff, default=to_json, encoding='utf-8', indent=2)
    finally:
        ff.close()

    main(['-d', '--basedir', basedir, '--datadir', 'data',
          '--run', 'config.txt'])
