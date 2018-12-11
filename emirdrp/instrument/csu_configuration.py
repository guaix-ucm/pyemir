#
# Copyright 2008-2018 Universidad Complutense de Madrid
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
from copy import deepcopy
import numpy as np

from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_MINIMUM_SLITLET_WIDTH_MM
from emirdrp.core import EMIR_MAXIMUM_SLITLET_WIDTH_MM


class CsuConfiguration(object):
    """Configurable Slit Unit (CSU) Configuration class definition.

    Attributes
    ----------
    _csu_bar_left : list of floats
        Location (mm) of the left bar for each slitlet.
    _csu_bar_right : list of floats
        Location (mm) of the right bar for each slitlet, using the
        same origin employed for _csu_bar_left (which is not the
        value stored in the FITS keywords.
    _csu_bar_slit_center : list of floats
        Middle point (mm) in between the two bars defining a slitlet.
    _csu_bar_slit_width : list of floats
        Slitlet width (mm), computed as the distance between the two
        bars defining the slitlet.

    """

    def __init__(self):
        self._csu_bar_left = None
        self._csu_bar_right = None
        self._csu_bar_slit_center = None
        self._csu_bar_slit_width = None

    def __str__(self):
        output = "<CsuConfiguration instance>\n"
        for i in range(EMIR_NBARS):
            ibar = i + 1
            strdum = "- [BAR{0:2d}] left, right, center, width: ".format(ibar)
            output += strdum
            strdum = "{0:7.3f} {1:7.3f} {2:7.3f} {3:7.3f}\n".format(
                self._csu_bar_left[i], self._csu_bar_right[i],
                self._csu_bar_slit_center[i], self._csu_bar_slit_width[i]
            )
            output += strdum
        return output

    def __eq__(self, other):
        if isinstance(other, CsuConfiguration):
            result = \
                (self._csu_bar_left == other._csu_bar_left) and \
                (self._csu_bar_right == other._csu_bar_right) and \
                (self._csu_bar_slit_center == other._csu_bar_slit_center) and \
                (self._csu_bar_slit_width == other._csu_bar_slit_width)
            return result
        return NotImplemented

    @classmethod
    def define_from_fits(cls, fitsobj, extnum=0):
        """Define class members from header information in FITS file.

        Parameters
        ----------
        fitsobj: file object
            FITS file whose header contains the CSU bar information
            needed to initialise the members of this class.
        extnum : int
            Extension number (first extension is 0)

        """

        # read input FITS file
        with fits.open(fitsobj) as hdulist:
            image_header = hdulist[extnum].header
            return cls.define_from_header(image_header)

    @classmethod
    def define_from_header(cls, image_header):
        """Define class members directly from FITS header.

        Parameters
        ----------
        image_header : instance of hdulist.header
            Header content from a FITS file.

        """

        self = CsuConfiguration()

        # declare lists to store configuration of CSU bars
        self._csu_bar_left = []
        self._csu_bar_right = []
        self._csu_bar_slit_center = []
        self._csu_bar_slit_width = []

        for i in range(EMIR_NBARS):
            ibar = i + 1
            keyword = 'CSUP{}'.format(ibar)
            if keyword in image_header:
                self._csu_bar_left.append(image_header[keyword])
            else:
                raise ValueError("Expected keyword " + keyword + " not found!")
            keyword = 'CSUP{}'.format(ibar + EMIR_NBARS)
            if keyword in image_header:
                # set the same origin as the one employed for _csu_bar_left
                self._csu_bar_right.append(341.5 - image_header[keyword])
            else:
                raise ValueError("Expected keyword " + keyword + " not found!")
            self._csu_bar_slit_center.append(
                (self._csu_bar_left[i] + self._csu_bar_right[i]) / 2
            )
            self._csu_bar_slit_width.append(
                self._csu_bar_right[i] - self._csu_bar_left[i]
            )

        return self

    def csu_bar_left(self, islitlet):
        """Return csu_bar_left for requested slitlet number."""

        return self._csu_bar_left[islitlet - 1]

    def csu_bar_right(self, islitlet):
        """Return csu_bar_right for requested slitlet number."""

        return self._csu_bar_right[islitlet - 1]

    def csu_bar_slit_center(self, islitlet):
        """Return csu_bar_slit_center for requested slitlet number."""

        return self._csu_bar_slit_center[islitlet - 1]

    def csu_bar_slit_width(self, islitlet):
        """Return csu_bar_slit_width for requested slitlet number."""

        return self._csu_bar_slit_width[islitlet - 1]

    def outdict(self, ndigits=3):
        """Return dictionary structure rounded to a given precision."""

        outdict = {}
        for i in range(EMIR_NBARS):
            ibar = i + 1
            cbar = 'slitlet' + str(ibar).zfill(2)
            outdict[cbar] = {}
            outdict[cbar]['_csu_bar_left'] = \
                round(self._csu_bar_left[i], ndigits)
            outdict[cbar]['_csu_bar_right'] = \
                round(self._csu_bar_right[i], ndigits)
            outdict[cbar]['_csu_bar_slit_center'] = \
                round(self._csu_bar_slit_center[i], ndigits)
            outdict[cbar]['_csu_bar_slit_width'] = \
                round(self._csu_bar_slit_width[i], ndigits)

        return outdict

    def widths_in_range_mm(
            self,
            minwidth=EMIR_MINIMUM_SLITLET_WIDTH_MM,
            maxwidth=EMIR_MAXIMUM_SLITLET_WIDTH_MM
    ):
        """Return list of slitlets which width is within given range

        Parameters
        ----------
        minwidth : float
            Minimum slit width (mm).
        maxwidth : float
            Maximum slit width (mm).

        Returns
        -------
        list_ok : list
            List of booleans indicating whether the corresponding
            slitlet width is within range

        """

        list_ok = []
        for i in range(EMIR_NBARS):
            slitlet_ok = minwidth <= self._csu_bar_slit_width[i] <= maxwidth
            if slitlet_ok:
                list_ok.append(i + 1)

        return list_ok


def merge_odd_even_csu_configurations(conf_odd, conf_even):
    """Merge CSU configuration using odd- and even-numbered values.

    The CSU returned CSU configuration include the odd-numbered values
    from 'conf_odd' and the even-numbered values from 'conf_even'.

    Parameters
    ----------
    conf_odd : CsuConfiguration instance
        CSU configuration corresponding to odd-numbered slitlets.
    conf_even : CsuConfiguration instance
        CSU configuration corresponding to even-numbered slitlets.

    Returns
    -------
    merged_conf : CsuConfiguration instance
        CSU configuration resulting from the merging process.

    """

    # initialize resulting CsuConfiguration instance using one of the
    # input configuration corresponding to the odd-numbered slitlets
    merged_conf = deepcopy(conf_odd)

    # update the resulting configuration with the values corresponding
    # to the even-numbered slitlets
    for i in range(EMIR_NBARS):
        ibar = i + 1
        if ibar % 2 == 0:
            merged_conf._csu_bar_left[i] = conf_even._csu_bar_left[i]
            merged_conf._csu_bar_right[i] = conf_even._csu_bar_right[i]
            merged_conf._csu_bar_slit_center[i] = \
                conf_even._csu_bar_slit_center[i]
            merged_conf._csu_bar_slit_width[i] = \
                conf_even._csu_bar_slit_width[i]

    # return merged configuration
    return merged_conf
