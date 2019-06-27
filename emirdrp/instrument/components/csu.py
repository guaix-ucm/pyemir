#
# Copyright 2019 Universidad Complutense de Madrid
#
# This file is part of EMIR DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from numina.instrument.hwdevice import HWDevice
from numina.array.bbox import BoundingBox
from numina.instrument.generic import ElementBase


class BarsModel(ElementBase):
    def __init__(self, name):
        super(BarsModel, self).__init__(name)


class ConfigurableSlitUnit(HWDevice):

    OPEN_L_LIMIT = 6
    OPEN_R_LIMIT = 2057
    NBARS = 55

    def __init__(self, name, barsmodel, parent=None):
        super(ConfigurableSlitUnit, self).__init__(name, parent=parent)

        # Attributtes added to read from JSON
        self.origin = None
        self.configurations = {}
        # Up to HERE

        self.lbars = []
        self.rbars = []
        self.slits = {}

        nominal = barsmodel.values['nominal']

        for i in range(self.NBARS):
            name = "BarL_{}".format(i + 1)
            _, ypos, delt, x0 = nominal[i]
            model = pol.Polynomial([x0, delt])
            bar = CSUBarL(name, ypos, model, parent=self)
            self.lbars.append(bar)

        for i in range(self.NBARS):
            name = "BarR_{}".format(self.NBARS + i + 1)
            _, ypos, delt, x0 = nominal[i]
            model = pol.Polynomial([x0, delt])
            bar = CSUBarR(name, ypos, model, parent=self)
            self.rbars.append(bar)

    # Methods added to read from JSON
    @classmethod
    def init_args(self, name, setup_objects):
        barsmodel = setup_objects['barsmodel']
        return (name, barsmodel), {}

    def is_open(self):
        lopen = all(bar.position <= self.OPEN_L_LIMIT for bar in self.lbars)
        ropen = all(bar.position >= self.OPEN_R_LIMIT for bar in self.rbars)

        return lopen and ropen

    def configure_with_header(self, hdr):

        for i in range(self.NBARS):
            csupos = hdr["CSUP{}".format(i + 1)]
            self.lbars[i].csupos = csupos
        for i in range(self.NBARS):
            csupos = hdr["CSUP{}".format(self.NBARS + i + 1)]
            self.lbars[i].csupos = csupos

    @property
    def test(self):
        return 0


class CSUBar(HWDevice):
    def __init__(self, name, ypos, model, parent=None, active=True):
        super(CSUBar, self).__init__(name, parent=parent)
        self.ypos = ypos
        self.csupos_ = 0
        self.active_ = active
        self.model = model
        # Geometry
        self.slit_h_virt = 16.242
        self.y0 = ypos
        self.y1 = ypos - self.slit_h_virt - 3
        self.y2 = ypos + self.slit_h_virt + 3
        self._x1 = 0
        self._x2 = 2047

    @property
    def csupos(self):
        return self.csupos_

    @csupos.setter
    def csupos(self, value):
        if self.active_:
            self.csupos_ = value
        else:
            # warning, bar is inactive
            pass

    @property
    def active(self):
        return self.active_

    @active.setter
    def active(self, value):
        self.active_ = value

    @property
    def position(self):
        return self.model(self.csupos_)

    @property
    def x1(self):
        return self._x1

    @property
    def x2(self):
        return self._x2

    def bbox(self):
        # origin is 1
        x1 = max(self.x1 - 1, 0)
        y1 = max(self.y1 - 1, 0)
        x2 = max(self.x2 - 1, 0, x1)
        y2 = max(self.y2 - 1, 0, y1)
        return BoundingBox.from_coordinates(x1, x2, y1, y2)


class CSUBarL(CSUBar):
    @property
    def x2(self):
        return self.position


class CSUBarR(CSUBar):
    @property
    def x1(self):
        return self.position


def create_test_header0():
    hdr = {}
    for i in range(55):
        hdr["CSUP{}".format(i + 1)] = 123
    for i in range(55, 110):
        hdr["CSUP{}".format(i + 1)] = 125
    return hdr


if __name__ == '__main__':

    import pkgutil
    import numpy
    import io
    import numpy.polynomial.polynomial as pol

    dumdata = pkgutil.get_data('emirdrp.instrument.configs', 'bars_nominal_positions_test.txt')
    ss = io.StringIO(dumdata.decode('utf8'))
    bars_nominal_positions = numpy.loadtxt(ss)

    class BarModel:
        def __init__(self):
            self.values = {}

    barmodel = BarModel()
    barmodel.values['nominal'] = bars_nominal_positions

    csu = ConfigurableSlitUnit("CSU", barmodel)

    print(csu.config_info())
    csu.configure({'CSU.BarL_1': {'csupos': 100}})
    print(csu.children)
    print(csu.is_open())
    dev = csu.get_device('CSU.BarL_12')
    csu.get_property('CSU.BarL_12.csupos')
    dev.csupos = 234
    print(dev.position)
    print(dev.bbox().slice)

    hdr = create_test_header0()

    csu.configure_with_header(hdr)
    print(csu.config_info())
