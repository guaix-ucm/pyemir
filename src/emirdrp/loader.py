#
# Copyright 2014-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import numina.core


def load_drp():
    return numina.core.drp_load('emirdrp', 'drp.yaml')
