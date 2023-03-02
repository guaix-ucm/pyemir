#
# Copyright 2008-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""The EMIR Data Reduction Pipeline"""

import logging

__version__ = '0.19'


# Top level NullHandler
logging.getLogger("emirdrp").addHandler(logging.NullHandler())
