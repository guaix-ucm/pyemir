#
# Copyright 2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from numina.store.load import load
from numina.core import Requirement
from numina.core import Result

import emirdrp.products as prods
from emirdrp.core.recipe import EmirRecipe


class RecWaveRecipe(EmirRecipe):
    """Builds a MasterRecWave from a seralized version"""
    filename = Requirement(str, 'Full path of MasterRecWave')
    master_rectwv = Result(prods.MasterRectWave)

    def run(self, rinput):
        filename = rinput.filename
        self.logger.debug('filename is %s', filename)
        master_rectwv = load(prods.MasterRectWave, filename)
        return self.create_result(
            master_rectwv=master_rectwv
        )
