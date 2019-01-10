#
# Copyright 2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from numina.core.insconf import ComponentFacade


class CsuComponent(ComponentFacade):
    def __init__(self, name, uuid, date_start, date_end=None, description=""):
        super(CsuComponent, self).__init__(name, uuid, date_start, date_end=date_end,
                                           description=description)
        self.configurations = {}

    def __setstate__(self, state):
        super(CsuComponent, self).__setstate__(state)
