#
# Copyright 2008-2012 Universidad Complutense de Madrid
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

import yaml
from numina.core import BaseInstrument, BasePipeline, InstrumentConfiguration
from numina.core import import_object
import logging
import pkgutil

import emir.recipes as recp
from emir.instrument.channels import convert_name_to_channels
from emir import __version__

_logger = logging.getLogger('emir')

def _load_modes():
    equiv_class = {}
    modes = [i for i in yaml.load_all(pkgutil.get_data('emir.instrument',
                                                   'obsmodes.yaml'))
          if i.instrument == 'EMIR']
    for m in modes:
        m.recipe_class = import_object(m.recipe)
        equiv_class[m.key] = m.recipe_class
    return modes, equiv_class

class EMIR_Pipeline(BasePipeline):
    def __init__(self, name, version, recipes):
        super(EMIR_Pipeline, self).__init__(name=name,
                    version=version,
                    recipes=recipes)

def sup1():

    _modes, equiv_class = _load_modes()

    _conf1 = InstrumentConfiguration('default', yaml.load(pkgutil.get_data('emir.instrument', 'default.yaml')))
    # converting
    _conf1 = convert_name_to_channels(_conf1)
    
    _conf2 = InstrumentConfiguration('alternate', yaml.load(pkgutil.get_data('emir.instrument', 'alt.yaml')))
    _conf2 = convert_name_to_channels(_conf2)

    class EMIR_Instrument(BaseInstrument):
        name = 'EMIR'
        configurations = {'default': _conf1,
                      'alt': _conf2}
        if _conf2 is None:
            del configurations['alt']
    
        pipelines = {'default': EMIR_Pipeline(name='default', version=__version__, recipes=equiv_class),
                     'alternate': EMIR_Pipeline(name='alternate', version=__version__, recipes=equiv_class)}
        modes = _modes

        def __init__(self, configuration):
            super(EMIR_Instrument, self).__init__('EMIR',
                    configuration)

    return EMIR_Instrument

def sup2():
    _conf = InstrumentConfiguration('default', yaml.load(pkgutil.get_data('emir.instrument', 'default.yaml')))

    class EMIR_MUX_Instrument(BaseInstrument):
        name = 'EMIR_MUX'
        configurations = {'default': _conf}
    
        pipelines = {'default': EMIR_Pipeline(name='default', version=__version__, recipes={})}
        modes = [i for i in yaml.load_all(pkgutil.get_data('emir.instrument',
                                                   'obsmodesmux.yaml'))
            if i.instrument == 'EMIR_MUX']

        def __init__(self, configuration):
            super(EMIR_MUX_Instrument, self).__init__('EMIR_MUX',
                    configuration)

    return EMIR_MUX_Instrument

