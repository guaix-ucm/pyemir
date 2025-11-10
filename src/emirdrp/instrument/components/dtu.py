#
# Copyright 2019-2025 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#
import astropy.wcs
import numpy

from numina.instrument.hwdevice import HWDevice


from .dtuaxis import DTUAxis
from .dtuaxis import apply_on_axis as apply_on_axis


def convert_md5(string):
    import hashlib

    mm = hashlib.md5()
    mm.update(string.encode("utf-8"))
    return mm.hexdigest()


class DetectorTranslationUnit(HWDevice):
    def __init__(self, name, origin=None, parent=None):
        super().__init__(name, origin=origin, parent=parent)

        # Attributes added to read from JSON
        self.origin = None
        self.configurations = {}
        # Up to HERE

        self.xaxis = DTUAxis("X", parent=self)
        self.yaxis = DTUAxis("Y", parent=self)
        self.zaxis = DTUAxis("Z", parent=self)

    @classmethod
    def from_component(
        cls,
        name: str,
        comp_id: str,
        origin=None,
        parent=None,
        properties=None,
        setup=None,
    ):
        obj = cls.__new__(cls)
        obj.__init__(comp_id, origin=origin, parent=parent)
        return obj

    def configure_with_header(self, hdr):
        self.xaxis.configure_with_header(hdr)
        self.yaxis.configure_with_header(hdr)
        self.zaxis.configure_with_header(hdr)
        self.is_configured = True

    def configure_with_img(self, img):
        """Create a DtuConf object from a FITS image"""
        if "MECS" in img:
            hdr = img["MECS"].header
        else:
            hdr = img[0].header
        self.configure_with_header(hdr)

    # everything defined with @property can be accessed
    # with get_property('name')

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
        dtur = self.coor_r
        xfac = dtur[0]
        yfac = -dtur[1]

        dtuwcs = astropy.wcs.WCS(naxis=2)
        dtuwcs.wcs.name = "DTU WCS"
        dtuwcs.wcs.crpix = [0, 0]
        dtuwcs.wcs.cdelt = [1, 1]
        dtuwcs.wcs.crval = [yfac, xfac]
        dtuwcs.wcs.ctype = ["linear", "linear"]
        dtuwcs.wcs.cunit = ["um", "um"]

        return dtuwcs

    def vector_shift(self, pixsize=None):
        import astropy.units as u
        import emirdrp.instrument.constants as cons

        if pixsize is None:
            pixsize = cons.EMIR_PIXSIZE
        # coordinates transformation from DTU coordinates
        # to image coordinates
        # Y inverted
        # XY switched
        # trans1 = [[1, 0, 0], [0,-1, 0], [0,0,1]]
        # trans2 = [[0,1,0], [1,0,0], [0,0,1]]
        trans3 = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]  # T3 = T2 * T1
        vec = numpy.dot(trans3, self.coor_r) * u.micron / pixsize
        return vec

    @classmethod
    def from_header(cls, hdr, name="DTU", parent=None, **kwds):
        """Create a DtuConf object from a FITS image"""
        obj = cls.__new__(cls)
        obj.__init__(name, parent=parent, **kwds)
        obj.configure_with_header(hdr)
        return obj

    @classmethod
    def from_img(cls, img, name="DTU", parent=None, **kwds):
        """Create a DtuConf object from a FITS image"""
        obj = cls.__new__(cls)
        obj.__init__(name, parent=parent, **kwds)
        obj.configure_with_img(img)
        return obj

    def allclose(self, other, rtol=1e-05, atol=1e-08, equal_nan=False):
        return (
            self.xaxis.allclose(other.xaxis, rtol=rtol, atol=atol, equal_nan=equal_nan)
            and self.yaxis.allclose(
                other.yaxis, rtol=rtol, atol=atol, equal_nan=equal_nan
            )
            and self.zaxis.allclose(
                other.zaxis, rtol=rtol, atol=atol, equal_nan=equal_nan
            )
        )

    def closeto(self, other, abserror):
        return self.allclose(other, rtol=0.0, atol=abserror, equal_nan=False)

    def stringify(self, ndig=3):
        xv = self.xaxis.stringify(ndig)
        yv = self.yaxis.stringify(ndig)
        zv = self.zaxis.stringify(ndig)
        return f"x={xv}:y={yv}:z={zv}"

    def __eq__(self, other):
        if isinstance(other, DtuConf):
            return (
                self.xaxis == other.xaxis
                and self.yaxis == other.yaxis
                and self.zaxis == other.zaxis
            )
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    def describe(self, ndig=3):
        # note: set the number of decimal figures in output to the precision
        # value employed in __eq__() function below

        output = (
            "<DtuConf instance>\n"
            "- XDTU..: {0.xaxis.coor:{width}.{prec}f}\n"
            "- YDTU..: {0.yaxis.coor:{width}.{prec}f}\n"
            "- ZDTU..: {0.zaxis.coor:{width}.{prec}f}\n"
            "- XDTU_0: {0.xaxis.coor_0:{width}.{prec}f}\n"
            "- YDTU_0: {0.yaxis.coor_0:{width}.{prec}f}\n"
            "- ZDTU_0: {0.zaxis.coor_0:{width}.{prec}f}\n"
            "- XDTU_F: {0.xaxis.coor_f:{width}.{prec}f}\n"
            "- YDTU_F: {0.yaxis.coor_f:{width}.{prec}f}\n"
            "- ZDTU_F: {0.zaxis.coor_f:{width}.{prec}f}\n"
        )
        return output.format(self, width=8, prec=ndig)

    def outdict(self, ndigits=3):
        """Return dictionary structure rounded to a given precision."""
        output = {
            "dtuhash": convert_md5(self.stringify(ndig=3)),
            "xdtu": round(self.xaxis.coor, ndigits),
            "xdtu_0": round(self.xaxis.coor_0, ndigits),
            "ydtu": round(self.yaxis.coor, ndigits),
            "ydtu_0": round(self.yaxis.coor_0, ndigits),
            "zdtu": round(self.zaxis.coor, ndigits),
            "zdtu_0": round(self.zaxis.coor_0, ndigits),
            "xdtu_f": round(self.xaxis.coor_f, ndigits),
            "ydtu_f": round(self.yaxis.coor_f, ndigits),
            "zdtu_f": round(self.zaxis.coor_f, ndigits),
        }
        return output


class DtuConf(DetectorTranslationUnit):
    pass


def average(*args):
    """Return instance with averaged values."""
    new = DtuConf("DTU_avg")
    for axis in ["xaxis", "yaxis", "zaxis"]:
        fun_res = apply_on_axis(numpy.mean, [getattr(arg, axis) for arg in args])
        getattr(new, axis).configure_me(fun_res)
    return new


def maxdiff(*args):
    """Return DtuConf instance with maximum differences."""
    new = DtuConf()
    for axis in ["xaxis", "yaxis", "zaxis"]:
        fun_res = apply_on_axis(numpy.ptp, [getattr(arg, axis) for arg in args])
        getattr(new, axis).configure_me(fun_res)
    return new
