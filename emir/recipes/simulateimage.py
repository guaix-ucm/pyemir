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

'''Recipe for image simulation.'''


import logging

import numpy

from numina.recipes import RecipeBase
import numina.qa as qa
from numina.simulation import run_counter
from numina.recipes.registry import ProxyQuery
from numina.recipes.registry import Schema
from emir.instrument.detector import Hawaii2Detector
from emir.dataproducts import create_raw
from emir.recipes import EmirRecipeMixin

_logger = logging.getLogger("emir.recipes")

class Recipe(RecipeBase, EmirRecipeMixin):
    '''Recipe to simulate EMIR images.'''
    
    required_parameters = [
        Schema('detector', ProxyQuery(dummy={}), 'Detector parameters'),
        Schema('readout', ProxyQuery(dummy={}), 'Readout mode'),
        Schema('name_format', 'r%05d', 'Filename format'),
    ]

    capabilities = ['simulate_image']
    
    def __init__(self, parameters, runinfo):
        super(Recipe, self).__init__(parameters, runinfo)
        #
        _logger.info('Run counter created')
        self.runcounter = run_counter("r%05d", suffix='')
        
        _logger.info('Creating detector')
        
        self.detector = Hawaii2Detector(**self.parameters['detector'])
        _logger.info('Configuring detector')
        self.detector.configure(self.parameters['readout'])
        
        self.input = numpy.zeros(self.parameters['detector']['shape'])
        self.detector.exposure(self.parameters['readout']['exposure'])
        self.repeat = self.parameters['readout']['repeat']
     
    def run(self):
        _logger.info('Creating simulated array')    
        output = self.detector.lpath(self.input)
        run = self.runcounter.next()
        headers = {'RUN': run, 'FILENAME': '%s.fits' % run}
        
        _logger.info('Collecting detector metadata')
        headers.update(self.detector.metadata())
        
        _logger.info('Building FITS structure')
        hdulist = create_raw(output, headers)
        return {'qa': qa.UNKNOWN, 'simulated_image': hdulist}
