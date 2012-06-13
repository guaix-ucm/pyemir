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

'''The EMIR Data Reduction Pipeline'''

import logging
import pkgutil
import importlib
    
import yaml
from numina.core import BaseInstrument, BasePipeline

import emir.recipes as recp
from emir import __version__


def super_import(path):
    spl = path.split('.')
    cls = spl[-1]
    mods = '.'.join(spl[:-1])
    mm = importlib.import_module(mods)
    Cls = getattr(mm, cls)
    return Cls

_modes = [i for i in yaml.load_all(pkgutil.get_data('emir.instrument',
                                                   'obsmodes.yaml'))
          if i.instrument == 'EMIR']

_equiv_class = {}

for _m in _modes:
    _m.recipe_class = super_import(_m.recipe)
    _equiv_class[_m.key] = _m.recipe_class

del _m

class EmirPipeline(BasePipeline):
    def __init__(self):
        super(EmirPipeline, self).__init__(name='EMIR', 
                version=__version__,
                recipes=_equiv_class)

class EMIR_Instrument(BaseInstrument):
    name = 'EMIR'
    modes = _modes
