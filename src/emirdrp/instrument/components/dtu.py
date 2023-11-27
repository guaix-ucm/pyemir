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
    """Represents one DTU axis of movement"""
    def __init__(self, name, parent=None, active=True):
        super().__init__(name, parent=parent)
        self.active_ = active
        self.coor_ = 0.0
        self.coor_f_ = 1.0
        self.coor_0_ = 0.0

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
    def coor_f(self):
        return self.coor_f_

    @coor_f.setter
    def coor_f(self, value):
        if self.active_:
            self.coor_f_ = value
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
    def coor_r(self):
        return (self.coor / self.coor_f) - self.coor_0

    def configure_with_header(self, hdr):
        self.coor = hdr[f'{self.name}DTU']
        self.coor_f = hdr.get('{}DTU_F'.format(self.name), 1.0)
        self.coor_0 = hdr.get('{}DTU_0'.format(self.name), 0.0)


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
