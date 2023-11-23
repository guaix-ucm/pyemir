#
# Copyright 2019-2023 Universidad Complutense de Madrid
#
# This file is part of EMIR DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from numina.instrument.hwdevice import HWDevice

class DTUAxis(HWDevice):
    def __init__(self, name, parent=None, active=True):
        super().__init__(name, parent=parent)
        self.active_ = active
        self.coor_ = 0
        self.coor_f_ = 0
        self.coor_0_ = 0

    @property
    def coor(self):
        return self.coor_

    @coor.setter
    def coor(self, value):
        if self.active_:
            self.coor_ = value
        else:
            # warning, bar is inactive
            pass

    @property
    def coor_0(self):
        return self.coor_0_

    @coor_0.setter
    def coor_0(self, value):
        if self.active_:
            self.coor_0_ = value
        else:
            # warning, bar is inactive
            pass

    @property
    def active(self):
        return self.active_

    @active.setter
    def active(self, value):
        self.active_ = value


class DetectorTranslationUnit(HWDevice):
    def __init__(self, name, parent=None, **kwds):
        super().__init__(name, parent=parent)

        # Attributtes added to read from JSON
        self.origin = None
        self.configurations = {}
        # Up to HERE

        self.xaxis = DTUAxis("X", parent=self)
        self.yaxis = DTUAxis("Y", parent=self)
        self.zaxis = DTUAxis("Z", parent=self)

    def configure_with_header(self, hdr):
        pass
