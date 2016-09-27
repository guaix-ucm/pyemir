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

"""Decorator for EMIR Recipes"""

from __future__ import print_function

import datetime
from numina.core import DataFrame, ObservationResult
import emirdrp.processing.info as info

def timeit(method):
    """Measure the time used by the recipe"""
    def timed_method(self, rinput):

        time_start = datetime.datetime.utcnow()
        result = method(self, rinput)
        time_end = datetime.datetime.utcnow()
        result.time_it(time_start, time_end)
        self.logger.info('total time measured')
        return result

    return timed_method


def loginfo(method):
    """Log the contents of Recipe Input"""


    def loginfo_method(self, rinput):

        klass = rinput.__class__

        for key in klass.stored():
            val = getattr(rinput, key)
            if isinstance(val, DataFrame):
                self.logger.debug("DataFrame %s", info.gather_info_dframe(val))
            elif isinstance(val, ObservationResult):
                for f in val.images:
                    self.logger.debug("OB DataFrame %s" , info.gather_info_dframe(f))
            else:
                pass

        result = method(self, rinput)

        return result

    return loginfo_method


def aggregate(method):

    def agregate_method(self, rinput):

        result = method(self, rinput)

        aggregate_m = getattr(self, 'aggregate_result')

        new_result = aggregate_m(self, result, rinput)

        return new_result

    return agregate_method
