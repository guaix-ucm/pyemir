#
# Copyright 2008-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Decorator for EMIR Recipes"""

from __future__ import print_function

from numina.core import DataFrame, ObservationResult
import emirdrp.processing.info as info


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
                    self.logger.debug("OB DataFrame %s",
                                      info.gather_info_dframe(f))
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
