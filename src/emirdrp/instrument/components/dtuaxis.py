#
# Copyright 2019-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import warnings

from astropy.io.fits import Header
from numina.instrument.hwdevice import HWDevice


class DTUAxis(HWDevice):
    """Represents one DTU axis of movement"""
    def __init__(self, name, parent=None, active=True):
        super().__init__(name, parent=parent)
        self.active_ = active
        self.coor_ = 0.0
        self.coor_f_ = 1.0
        self.coor_0_ = 0.0

    def allclose(self, other, rtol=1e-05, atol=1e-08, equal_nan=False):
        import numpy
        a = [self.coor, self.coor_f, self.coor_0]
        b = [other.coor, other.coor_f, other.coor_0]
        return numpy.allclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)

    def closeto(self, other, abserror):
        return self.allclose(other, rtol=0.0, atol=abserror, equal_nan=False)

    def set_prop_key(self, key, value):
        if self.active_:
            setattr(self, key, value)
        else:
            warnings.warn("DTUAxis is inactive")

    @property
    def coor(self):
        return self.coor_

    @coor.setter
    def coor(self, value):
        self.set_prop_key('coor_', value)

    @property
    def coor_0(self):
        return self.coor_0_

    @coor_0.setter
    def coor_0(self, value):
        self.set_prop_key('coor_0_', value)

    @property
    def coor_f(self):
        return self.coor_f_

    @coor_f.setter
    def coor_f(self, value):
        self.set_prop_key('coor_f_', value)

    @property
    def active(self):
        return self.active_

    @active.setter
    def active(self, value):
        self.active_ = value

    @property
    def coor_r(self):
        return (self.coor / self.coor_f) - self.coor_0

    def configure_me_with_header(self, hdr):
        if isinstance(hdr, dict):
            hdr = Header(hdr)
        super().configure_me_with_header(hdr)
        self.coor = hdr[f'{self.name}DTU']
        self.coor_f = hdr.get(f'{self.name}DTU_F', 1.0)
        self.coor_0 = hdr.get(f'{self.name}DTU_0', 0.0)
        # FIXME: this is not needed
        self.is_configured = True

    def configure_me_with_image(self, image):
        if 'MECS' in image:
            hdr = image['MECS'].header
        else:
            hdr = image[0].header
        self.configure_me_with_header(hdr)

    @classmethod
    def from_header(cls, hdr, name, parent=None):
        obj = cls.__new__(cls)
        obj.__init__(name, parent=parent)
        obj.configure_with_header(hdr)
        return obj

    def stringify(self, ndig=3):
        values = (self.coor_, self.coor_f_, self.coor_0_)
        m1 = [round(r, ndig) for r in values]
        return f'{m1[0]:+010.3f}:{m1[1]:+010.3f}:{m1[2]:+010.3f}'

    def __hash__(self):
        return hash(self.stringify())

    def __eq__(self, other):
        ndig = 3
        if isinstance(other, DTUAxis):
            result = \
                self.name.lower() == other.name.lower() and \
                self.stringify(ndig) == other.stringify(ndig)

            return result
        return NotImplemented

    def __ne__(self, other):
        return not self == other


def apply_on_axis(func, args):
    """Apply a function to the members of a sequence of DtuAxis"""
    coor_avg = float(func([arg.coor for arg in args]))
    coor_f_avg = float(func([arg.coor_f for arg in args]))
    coor_0_avg = float(func([arg.coor_0 for arg in args]))

    return dict(coor=coor_avg, coor_f=coor_f_avg, coor_0=coor_0_avg)
