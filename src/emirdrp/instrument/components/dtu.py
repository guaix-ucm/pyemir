#
# Copyright 2019-2023 Universidad Complutense de Madrid
#
# This file is part of EMIR DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from numina.instrument.hwdevice import HWDevice


from .dtuaxis import DTUAxis

class DetectorTranslationUnit(HWDevice):
    def __init__(self, name, parent=None, **kwds):
        super().__init__(name, parent=parent)

        # Attributes added to read from JSON
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
