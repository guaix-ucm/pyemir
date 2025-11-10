#
# Copyright 2008-2024 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Decorator for EMIR Recipes"""

from numina.core import DataFrame, ObservationResult
from emirdrp.datamodel import EmirDataModel


def loginfo(method):
    """Log the contents of Recipe Input"""

    # FIXME: datamodel is probably already defined in method
    # or in rinput
    datamodel = EmirDataModel()

    def loginfo_method(self, rinput):

        klass = rinput.__class__
        for key in klass.stored():
            val = getattr(rinput, key)
            if isinstance(val, DataFrame):
                self.logger.debug("DataFrame %s", datamodel.gather_info_dframe(val))
            elif isinstance(val, ObservationResult):
                for msg in datamodel.gather_info_oresult(val):
                    self.logger.debug("OB DataFrame %s", msg)
            else:
                pass

        result = method(self, rinput)

        return result

    return loginfo_method


def aggregate(method):

    def agregate_method(self, rinput):

        result = method(self, rinput)

        aggregate_m = getattr(self, "aggregate_result")

        new_result = aggregate_m(self, result, rinput)

        return new_result

    return agregate_method
