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

'''Intensity Flatfield Recipe.'''

import logging
import os

import numpy
import pyfits

from numina.recipes import RecipeBase
from numina.image import DiskImage
from numina.image.flow import SerialFlow
from numina.image.processing import DarkCorrector, NonLinearityCorrector, BadPixelCorrector
from numina.array.combine import flatcombine
from numina.worker import para_map
from numina.recipes.registry import ProxyPath, ProxyQuery
from numina.recipes.registry import Schema
import numina.qa as QA
from emir.dataproducts import create_result
from emir.recipes import EmirRecipeMixin

_logger = logging.getLogger("emir.recipes")

class Recipe(RecipeBase, EmirRecipeMixin):
    '''Recipe to process data taken in intensity flat-field mode.
        
    Recipe to process intensity flat-fields. The flat-on and flat-off images are
    combined (method?) separately and the subtracted to obtain a thermal subtracted
    flat-field.
    
    **Observing modes:**
    
     * Intensity Flat-Field
    
    **Inputs:**
    
      * A list of lamp-on flats
      * A list of lamp-off flats
      * A master dark frame
      * A model of the detector. 
    
    **Outputs:**
    
     * TBD
    
    **Procedure:**
    
     * A combined thermal subtracted flat field, normalized to median 1, 
       with with variance extension and quality flag. 
    
    '''
    
    required_parameters = [
        Schema('images', ProxyPath('/observing_block/result/images'), 'A list of paths to images'),
        Schema('master_bias', ProxyQuery(), 'Master bias image'),
        Schema('master_dark', ProxyQuery(), 'Master dark image'),
        Schema('master_bpm', ProxyQuery(), 'Master bad pixel mask'),
        Schema('nonlinearity', ProxyQuery(dummy=[1.0, 0.0]), 'Polynomial for non-linearity correction'),
        Schema('output_filename', 'result.fits', 'Name of the output image')
    ]
    
    capabilities = ['intensity_flatfield']
    
    def __init__(self, param, runinfo):
        super(Recipe, self).__init__(param, runinfo)
                        
    def basic_processing(self, image, flow):
        hdulist = pyfits.open(image)
        try:
            hdu = hdulist['primary']                
            # Processing
            hdu = flow(hdu)
            
            hdulist.writeto(os.path.basename(image), clobber=True)
        finally:
            hdulist.close()
    
    def compute_superflat(self):
        try:
            filelist = []
            data = []
            for image in self.parameters['images']:
                hdulist = pyfits.open(image.filename, memmap=True, mode='readonly')
                filelist.append(hdulist)
                data.append(hdulist['primary'].data)
            _logger.info('Computing scale factors')            
            scales = [numpy.median(datum) for datum in data]
        
            _logger.info('Combining')
            illum = flatcombine(data, scales=scales, method='median', 
                                blank=1.0 / scales[0])
            return illum
        finally:
            for fileh in filelist:               
                fileh.close()        

    
    
    def setup(self):
        self.parameters['master_dark'] = DiskImage(os.path.abspath(self.parameters['master_dark']))
        self.parameters['master_bpm'] = DiskImage(os.path.abspath(self.parameters['master_bpm']))
            
        self.parameters['images'] = [DiskImage(os.path.abspath(path)) 
                                     for path in self.parameters['images']]            

    def run(self):
        #nthreads = self.runinfo['nthreads']
        primary_headers = {'FILENAME': self.parameters['output_filename']}

        bpm = pyfits.getdata(self.parameters['master_bpm'].filename)
        dark = pyfits.getdata(self.parameters['master_dark'].filename)
        
        basicflow = SerialFlow([BadPixelCorrector(bpm),
                           NonLinearityCorrector(self.parameters['nonlinearity']),
                           DarkCorrector(dark)])
        
        for image in self.parameters['images']:
            self.basic_processing(image.filename, basicflow)        
        
        # Illumination seems to be necessary
        # ----------------------------------
        
        illum_data = self.compute_superflat()
            
    
        _logger.info("Final image created")
        illum = create_result(illum_data[0], headers=primary_headers,
                                   variance=illum_data[1],
                                   exmap=illum_data[2])
        
        return {'qa': QA.UNKNOWN, 'illumination_image': illum}
    
