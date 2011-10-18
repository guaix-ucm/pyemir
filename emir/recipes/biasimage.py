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

'''Bias image recipe and associates.'''

import logging
import os

from numina.array.combine import zerocombine
from numina.image import DiskImage
import numina.qa
from numina.recipes import RecipeBase
from emir.dataproducts import create_result
from emir.recipes import EmirRecipeMixin
from numina.recipes.registry import ProxyPath
from numina.recipes.registry import Schema

_logger = logging.getLogger("emir.recipes")

class Recipe(RecipeBase, EmirRecipeMixin):
    '''Recipe to process data taken in Bias image Mode.

    Bias images only appear in Simple Readout mode.

    **Observing modes:**
    
     * Bias Image (3.1)   
    
    **Inputs:**
    
     * A list of bias images
     * A model of the detector (gain, RN)
    
    **Outputs:**
    
     * A combined bias frame, with variance extension and quality flag. 
    
    **Procedure:**
    
    The list of images can be readly processed by combining them with a typical
    sigma-clipping algorithm.
    
    '''
    
    required_parameters = [
        Schema('images', ProxyPath('/observing_block/result/images'), 'A list of paths to bias images'),
        Schema('combine', 'median', 'Combine method'),
        Schema('output_filename', 'bias.fits', 'Name of the bias output image')
    ]

    capabilities = ['bias_image']
    
    def __init__(self, param, runinfo):
        super(Recipe, self).__init__(param, runinfo)
        
    def  setup(self):
        # Sanity check, check: all images belong to the same detector mode
        self.parameters['images'] = [DiskImage(os.path.abspath(path)) 
                  for path in self.parameters['images']]
        
    def run(self):        
        primary_headers = {'FILENAME': self.parameters['output_filename']}
        
        images = self.parameters['images']
        
        # Open all zero images
        alldata = []
        allmasks = []
        try:
            for n in images:
                n.open(mode='readonly', memmap=True)
                alldata.append(n.data)
            
            # Combine them
            cube = zerocombine(alldata, allmasks, 
                               method=self.parameters['combine'])
        
            result = create_result(cube[0], headers=primary_headers,
                                   variance=cube[1], 
                                   exmap=cube[2].astype('int16'))
        
            return {'qa': numina.qa.UNKNOWN, 'bias_image': result}
        finally:
            for n in images:
                n.close()

if __name__ == '__main__':
    import os
    import uuid
    
    from numina.compatibility import json
    from numina.user import main
    from numina.jsonserializer import to_json

    logging.basicConfig(level=logging.DEBUG)
    _logger.setLevel(logging.DEBUG)
        
    pv = {'recipes': {'default': {'parameters': {
                                'output_filename': 'perryr.fits',
                                'combine': 'average',
                                },
                     'run': {'repeat': 1,
                             },   
                        }},
          'observing_block': {'instrument': 'emir',
                       'mode': 'bias_image',
                       'id': 1,
                       'result': {
                                  'images': ['apr21_0046.fits', 
                                             'apr21_0047.fits', 
                                             'apr21_0048.fits',
                                             'apr21_0049.fits', 
                                             'apr21_0050.fits'
                                             ],
                            },
                       },                     
    }
    
    os.chdir('/home/spr/Datos/emir/test7')
    
    
    # Creating base directory for storing results
    uuidstr = str(uuid.uuid1()) 
    basedir = os.path.abspath(uuidstr)
    os.mkdir(basedir)
    
    ff = open('config.txt', 'w+')
    try:
        json.dump(pv, ff, default=to_json, encoding='utf-8', indent=2)
    finally:
        ff.close()

    main(['-d','--basedir', basedir, '--datadir', 'data', '--run', 'config.txt'])

