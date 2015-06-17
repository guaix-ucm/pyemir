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

from numina.store import dump, load

from .dataproducts import ChannelLevelStatistics

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
