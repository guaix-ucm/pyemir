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

import yaml
from numina.pipeline import register_pipeline
from numina.pipeline import BaseInstrument

from emir.recipes import EmirPipeline
import emir.recipes as recp
from emir import __version__

register_pipeline(EmirPipeline(version=__version__))

_modes = [i for i in yaml.load_all(pkgutil.get_data('emir.instrument',
                                                   'obsmodes.yaml'))
          if i.instrument == 'EMIR']

for m in _modes:
    _fqn = m.recipe.split('.')[-1]
    m.recipe_class = getattr(recp, _fqn)

del m

class EMIR_Instrument(BaseInstrument):
    name = 'EMIR'
    modes = _modes
