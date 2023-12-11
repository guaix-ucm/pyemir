#
# Copyright 2008-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""The EMIR Data Reduction Pipeline"""

import logging

import emirdrp._version

__version__ = emirdrp._version.__version__


# Top level NullHandler
logging.getLogger("emirdrp").addHandler(logging.NullHandler())
