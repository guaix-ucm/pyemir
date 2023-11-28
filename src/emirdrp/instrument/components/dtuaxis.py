#
# Copyright 2019-2023 Universidad Complutense de Madrid
#
# This file is part of EMIR DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

import contextlib
import warnings

from numina.instrument.hwdevice import HWDevice


class DTUAxis(HWDevice):
    """Represents one DTU axis of movement"""
    def __init__(self, name, parent=None, active=True):
        super().__init__(name, parent=parent)
        self.active_ = active
        self.coor_ = 0.0
        self.coor_f_ = 1.0
        self.coor_0_ = 0.0

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
        super().configure_me_with_header(hdr)
        self.coor = hdr[f'{self.name}DTU']
        self.coor_f = hdr.get('{}DTU_F'.format(self.name), 1.0)
        self.coor_0 = hdr.get('{}DTU_0'.format(self.name), 0.0)

    def configure_me_with_imge(self, image):
        if 'MECS' in image:
            hdr = image['MECS'].header
        else:
            hdr = image[0].header
        self.configure_me_with_header(hdr)


@contextlib.contextmanager
def managed_ndig(obj, ndig):
    """Context manager to handle ndig"""
    old_ndig = obj.get_ndig()
    obj.set_ndig(ndig)
    try:
        yield ndig
    finally:
        # Code to release resource, e.g.:
        obj.set_ndig(old_ndig)


class DtuAxisAdaptor(DTUAxis):
    """Represents one DTU axis of movement"""

    def __init__(self, name, coor=0.0, coor_f=1.0, coor_0=0.0):
        super().__init__(name)
        if name.lower() not in ['x', 'y', 'z']:
            raise ValueError('"name" must be "X", "Y" or "Z')

        # Parent values
        self.coor = coor
        self.coor_f = coor_f
        self.coor_0 = coor_0
        # Upto here
        self.ndig = 3

    @classmethod
    def from_header(cls, hdr, name):
        if name.lower() not in ['x', 'y', 'z']:
            raise ValueError('"name" must be "X", "Y" or "Z')
        coor = hdr['{}DTU'.format(name)]
        coor_f = hdr.get('{}DTU_F'.format(name), 1.0)
        coor_0 = hdr.get('{}DTU_0'.format(name), 0.0)
        return DtuAxisAdaptor(name, coor, coor_f, coor_0)

    def allclose(self, other, rtol=1e-05, atol=1e-08, equal_nan=False):
        import numpy
        a = [self.coor, self.coor_f, self.coor_0]
        b = [other.coor, other.coor_f, other.coor_0]
        return numpy.allclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)

    def closeto(self, other, abserror):
        return self.allclose(other, rtol=0.0, atol=abserror, equal_nan=False)

    def set_ndig(self, ndig):
        self.ndig = ndig

    def get_ndig(self):
        return self.ndig

    def __eq__(self, other):
        # note: set the precision (number of decimal places) to the same
        # number employed in __str__() function above to print out member
        # values
        # FIXME: this is eqv to allclose with atol=1e-4?
        if isinstance(other, DtuAxisAdaptor):
            result = \
                self.name.lower() == other.name.lower() and \
                (round(self.coor, self.ndig) == round(other.coor, self.ndig)) and \
                (round(self.coor_f, self.ndig) == round(other.coor_f, self.ndig)) and \
                (round(self.coor_0, self.ndig) == round(other.coor_0, self.ndig))

            return result
        return NotImplemented

    def __ne__(self, other):
        return not self == other


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
        self.xaxis.configure_with_header(hdr)
        self.yaxis.configure_with_header(hdr)
        self.zaxis.configure_with_header(hdr)

    def configure_with_img(self, img):
        """Create a DtuConf object from a FITS image"""
        if 'MECS' in img:
            hdr = img['MECS'].header
        else:
            hdr = img[0].header
        self.configure_with_header(hdr)

    @property
    def xdtu_r(self):
        return self.xaxis.coor_r

    @property
    def ydtu_r(self):
        return self.yaxis.coor_r

    @property
    def zdtu_r(self):
        return self.zaxis.coor_r

    @property
    def coor_r(self):
        return [self.xaxis.coor_r, self.yaxis.coor_r, self.zaxis.coor_r]

    @property
    def coor(self):
        return [self.xaxis.coor, self.yaxis.coor, self.zaxis.coor]

    def to_wcs(self):
        """Create a WCS structure for DTU measurements"""
        import astropy.wcs
        dtur = self.coor_r
        xfac = dtur[0]
        yfac = -dtur[1]

        dtuwcs = astropy.wcs.WCS(naxis=2)
        dtuwcs.wcs.name = 'DTU WCS'
        dtuwcs.wcs.crpix = [0, 0]
        dtuwcs.wcs.cdelt = [1, 1]
        dtuwcs.wcs.crval = [yfac, xfac]
        dtuwcs.wcs.ctype = ['linear', 'linear']
        dtuwcs.wcs.cunit = ['um', 'um']

        return dtuwcs

    @classmethod
    def from_header(cls, hdr, name="DTU"):
        """Create a DtuConf object from a FITS image"""
        obj = DetectorTranslationUnit(name=name)
        obj.configure_with_header(hdr)
        return obj

    @classmethod
    def from_img(cls, img, name="DTU"):
        """Create a DtuConf object from a FITS image"""
        obj = DetectorTranslationUnit(name=name)
        obj.configure_with_img(img)
        return obj
