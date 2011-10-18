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

'''Dark image recipe.''' 

import os
import logging

from numina.recipes import RecipeBase, RecipeError
from numina.array.combine import mean
from numina.image import DiskImage
import numina.qa
from emir.dataproducts import create_result
from emir.recipes import EmirRecipeMixin
from numina.recipes.registry import ProxyPath
from numina.recipes.registry import Schema

_logger = logging.getLogger("emir.recipes")

class Recipe(RecipeBase, EmirRecipeMixin):
    '''Recipe to process data taken in Dark current image Mode.

    Recipe to process dark images. The dark images will be combined 
    weighting with the inverses of the corresponding variance extension. 
    They do not have to be of the same exposure time t, they will be 
    scaled to the same t0 ~ 60s (those with very short exposure time 
    should be avoided).
    
    **Observing mode:**
    
     * Dark current Image (3.2)
    
    **Inputs:**
    
     * A list of dark images 
     * A model of the detector (gain, RN)
    
    **Outputs:**
    
     * A combined dark frame, with variance extension and quality flag. 
    
    **Procedure:**
    
    The process will remove cosmic rays (using a typical sigma-clipping algorithm).
    
    ''' 

    capabilities = ['dark_current_image']
    
    required_parameters = [
        Schema('images', ProxyPath('/observing_block/result/images'), 'A list of paths of dark images'),
        Schema('output_filename', 'dark.fits', 'Name of the dark output image')
    ]    
        
    def __init__(self, param, runinfo):
        super(Recipe, self).__init__(param, runinfo)
        # Default parameters. This can be read from a file        
    
    def  setup(self):
        
        self.parameters['images'] = [DiskImage(os.path.abspath(path)) 
                  for path in self.parameters['images']]
        
        # Sanity check, check: all images belong to the same detector mode
        def readmode(hdr):
            rom = [hdr[key] for key in ['EXPOSED','READSCHM', 
                                         'READMODE', 'READNUM', 'READREPT', 
                                         'DARKTIME', 'EXPTIME']]
            return rom
        
        first = self.parameters['images'][0]
        first.open(mode='readonly')
        hdr = first.meta
        first_rom = readmode(hdr)
        first.close()
        for image in self.parameters['images'][1:]:
            image.open(mode='readonly')
            hdr = image.meta
            rom = readmode(hdr)
            _logger.debug('Readmode parameters from image %s are %s', image.filename, rom)
            if rom != first_rom:
                _logger.warning('Got %s, expected %s', rom, first_rom)
                raise RecipeError
            image.close()
           
    def run(self):
        primary_headers = {'FILENAME': self.parameters['output_filename']}
        
        images = self.parameters['images']

        alldata = []
        allmasks = []
        
        try:
            for n in images:
                n.open(mode='readonly', memmap=True)
                alldata.append(n.data)
            
            # Combine them
            cube = mean(alldata, allmasks)
        
            result = create_result(cube[0], headers=primary_headers, 
                                   variance=cube[1], 
                                   exmap=cube[2].astype('int16'))
        
            return {'qa': numina.qa.UNKNOWN, 'dark_image': result}

        finally:
            for n in images:
                n.close()

    