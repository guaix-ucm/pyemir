#
# Copyright 2008-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyEmir is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
#

from __future__ import division
from __future__ import print_function

from astropy.io import fits
import numpy as np


class DtuConfiguration(object):
    """Detector Translation Unit (DTU) Configuration class definition.

    Attributes
    ----------
    xdtu : float
        XDTU fits keyword value.
    ydtu : float
        YDTU fits keyword value.
    zdtu : float
        ZDTU fits keyword value.
    xdtu_0 : float
        XDTU_0 fits keyword value.
    ydtu_0 : float
        YDTU_0 fits keyword value.
    zdtu_0 : float
        ZDTU_0 fits keyword value.

    """

    def __init__(self):
        self.xdtu = None
        self.ydtu = None
        self.zdtu = None
        self.xdtu_0 = None
        self.ydtu_0 = None
        self.zdtu_0 = None

    def __str__(self):
        # note: set the number of decimal figures in output to the precision
        # value employed in __eq__() function below
        output = "<DtuConfiguration instance>\n"
        strdum = "- XDTU..: {0:8.3f}\n".format(self.xdtu)
        output += strdum
        strdum = "- YDTU..: {0:8.3f}\n".format(self.ydtu)
        output += strdum
        strdum = "- ZDTU..: {0:8.3f}\n".format(self.zdtu)
        output += strdum
        strdum = "- XDTU_0: {0:8.3f}\n".format(self.xdtu_0)
        output += strdum
        strdum = "- YDTU_0: {0:8.3f}\n".format(self.ydtu_0)
        output += strdum
        strdum = "- ZDTU_0: {0:8.3f}".format(self.zdtu_0)
        output += strdum
        return output

    def __eq__(self, other):
        # note: set the precision (number of decimal places) to the same
        # number employed in __str__() function above to print out member
        # values
        if isinstance(other, DtuConfiguration):
            ndig = 3
            result = \
                (round(self.xdtu, ndig) == round(other.xdtu, ndig)) and \
                (round(self.ydtu, ndig) == round(other.ydtu, ndig)) and \
                (round(self.zdtu, ndig) == round(other.zdtu, ndig)) and \
                (round(self.xdtu_0, ndig) == round(other.xdtu_0, ndig)) and \
                (round(self.ydtu_0, ndig) == round(other.ydtu_0, ndig)) and \
                (round(self.zdtu_0, ndig) == round(other.zdtu_0, ndig))
            return result
        return NotImplemented

    def __ne__(self, other):
        return not self == other

    @classmethod
    def define_from_fits(cls, fitsobj, extnum=0):
        """Define class object from header information in FITS file.

        Parameters
        ----------
        fitsobj: file object
            FITS file whose header contains the DTU information
            needed to initialise the members of this class.
        extnum : None or int or str
            Extension number or extension name. If None, select 'MECS' extension
            if available. If not, use primary extension

        """

        # read input FITS file
        with fits.open(fitsobj) as hdulist:
            if extnum is None:
                if 'MECS' in hdulist:
                    extnum = 'MECS'
                else:
                    extnum = 0
            image_header = hdulist[extnum].header
            return cls.define_from_header(image_header)

    @classmethod
    def define_from_img(cls, img):
        """Define class object from header information in HDUList object.

        Parameters
        ----------
        img: HDUList
            FITS image whose header contains the DTU information
            needed to initialise the members of this class.

        """

        # read input FITS file
        if 'MECS' in img:
            extnum = 'MECS'
        else:
            extnum = 0
        image_header = img[extnum].header
        return cls.define_from_header(image_header)

    @classmethod
    def define_from_header(cls, image_header):
        """Define class object directly from FITS header.

        Parameters
        ----------
        image_header : instance of fits.Header
            Header content from a FITS file.

        """
        img_h_u = dict((k.lower(), v) for k, v in image_header.items())
        return cls.define_from_dictionary(img_h_u)

    @classmethod
    def define_from_dictionary(cls, inputdict):
        """Define class object from dictionary.

        Parameters
        ----------
        inputdict : dictionary like object
            Dictionary like object defining each member of the class.

        """

        self = DtuConfiguration()
        for item in self.__dict__:
            self.__dict__[item] = inputdict[item]
        return self

    @classmethod
    def define_from_values(cls, xdtu, ydtu, zdtu, xdtu_0, ydtu_0, zdtu_0):
        """Define class object from from provided values.

        Parameters
        ----------
        xdtu : float
            XDTU fits keyword value.
        ydtu : float
            YDTU fits keyword value.
        zdtu : float
            ZDTU fits keyword value.
        xdtu_0 : float
            XDTU_0 fits keyword value.
        ydtu_0 : float
            YDTU_0 fits keyword value.
        zdtu_0 : float
            ZDTU_0 fits keyword value.

        """

        self = DtuConfiguration()
        # define DTU variables
        self.xdtu = xdtu
        self.ydtu = ydtu
        self.zdtu = zdtu
        self.xdtu_0 = xdtu_0
        self.ydtu_0 = ydtu_0
        self.zdtu_0 = zdtu_0
        return self

    def closeto(self, other, abserror):
        """Check that all the members are equal within provided absolute error.

        Parameters
        ----------
        other : DtuConfiguration object
            DTU configuration instance to be compared with self.
        abserror : float
            Absolute maximum allowed error.

        Returns
        -------
        result : bool
            True is all members are within the specified maximum
            absolute error

        """

        result = \
            (abs(self.xdtu - other.xdtu) <= abserror) and \
            (abs(self.ydtu - other.ydtu) <= abserror) and \
            (abs(self.zdtu - other.zdtu) <= abserror) and \
            (abs(self.xdtu_0 - other.xdtu_0) <= abserror) and \
            (abs(self.ydtu_0 - other.ydtu_0) <= abserror) and \
            (abs(self.zdtu_0 - other.zdtu_0) <= abserror)
        return result

    def outdict(self, ndigits=3):
        """Return dictionary structure rounded to a given precision."""

        output = self.__dict__.copy()
        for item in output:
            output[item] = round(output[item], ndigits)
        return output


def average_dtu_configurations(list_of_objects):
    """Return DtuConfiguration instance with averaged values.

    Parameters
    ----------
    list_of_objects : python list
        List of DtuConfiguration instances to be averaged.

    Returns
    -------
    result : DtuConfiguration instance
        Object with averaged values.

    """

    result = DtuConfiguration()

    if len(list_of_objects) == 0:
        return result

    list_of_members = result.__dict__.keys()

    # compute average of all the members of the class
    for member in list_of_members:
        result.__dict__[member] = np.mean(
            [tmp_dtu.__dict__[member] for tmp_dtu in list_of_objects]
        )

    return result


def maxdiff_dtu_configurations(list_of_objects):
    """Return DtuConfiguration instance with maximum differences.

    Parameters
    ----------
    list_of_objects : python list
        List of DtuConfiguration instances to be averaged.

    Returns
    -------
    result : DtuConfiguration instance
        Object with averaged values.

    """

    result = DtuConfiguration()

    if len(list_of_objects) == 0:
        return result

    list_of_members = result.__dict__.keys()

    # compute maximum difference for each member
    for member in list_of_members:
        tmp_array = np.array(
            [tmp_dtu.__dict__[member] for tmp_dtu in list_of_objects]
        )
        minval = tmp_array.min()
        maxval = tmp_array.max()
        result.__dict__[member] = maxval - minval

    return result
