#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

from numina.user.cli import main


def run_recipe():
    main(['run', 'obsrun.yaml', '-r', 'control.yaml'])
