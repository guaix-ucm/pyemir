#
# Copyright 2008-2014 Universidad Complutense de Madrid
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

'''Data products produced by the EMIR pipeline.'''

import logging

import numpy
import yaml

from numina.store import dump, load

from .products import ChannelLevelStatistics
from .products import LinesCatalog
from .products import SlitsCatalog
from emirdrp.wavecal.slitlet import Slitlet

_logger = logging.getLogger('emirdrp.store')

_logger.debug('register dump functions')

@dump.register(ChannelLevelStatistics)
def _(tag, obj, where):
    fname = 'statistics.txt'

    header = '''Channel Level Statistics
comment 2
pixels start in 1
pixels end in 2048
exposure=%s
xbegin xend ybegin yend mean median var

'''
    inter = header % obj.exposure
    numpy.savetxt(fname, obj.statistics, header=inter)
    return fname

_logger.debug('register load functions')

@load.register(LinesCatalog)
def _l(tag, obj):

    with open(obj, 'r') as fd:
        linecat = numpy.genfromtxt(fd)
    return linecat


@load.register(SlitsCatalog)
def _l(tag, obj):

    with open(obj, 'r') as fd:
        slits_cat = yaml.load(fd)

    slits_list = []
    for slit in slits_cat:
        bbox = slit['bbox']
        slitdum = Slitlet(*bbox)
        borders = slit['borders']
        slitdum.set_nc_coeff_lower_boundary_pix(borders[0])
        slitdum.set_nc_coeff_upper_boundary_pix(borders[1])
        slits_list.append(slitdum)
    return slits_list
