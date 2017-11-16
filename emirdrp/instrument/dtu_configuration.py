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
    defined : bool
        Indicates whether the DTU parameters have been properly defined.

    """

    def __init__(self):
        self.xdtu = None
        self.ydtu = None
        self.zdtu = None
        self.xdtu_0 = None
        self.ydtu_0 = None
        self.zdtu_0 = None
        self.defined = False

    def __str__(self):
        # note: set the number of decimal figures in output to the precision
        # value employed in __eq__() function below
        output = "<DtuConfiguration instance>\n"
        if self.defined:
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
        else:
            output += "- XDTU..:  None\n"
            output += "- YDTU..:  None\n"
            output += "- ZDTU..:  None\n"
            output += "- XDTU_0:  None\n"
            output += "- YDTU_0:  None\n"
            output += "- ZDTU_0:  None"
        return output

    def __eq__(self, other):
        # note: set the precision (number of decimal places) to the same
        # number employed in __str__() function above to print out member
        # values
        ndig = 3
        result = \
            (self.defined == other.defined) and \
            (round(self.xdtu, ndig) == round(other.xdtu, ndig)) and \
            (round(self.ydtu, ndig) == round(other.ydtu, ndig)) and \
            (round(self.zdtu, ndig) == round(other.zdtu, ndig)) and \
            (round(self.xdtu_0, ndig) == round(other.xdtu_0, ndig)) and \
            (round(self.ydtu_0, ndig) == round(other.ydtu_0, ndig)) and \
            (round(self.zdtu_0, ndig) == round(other.zdtu_0, ndig))
        return result

    def __ne__(self, other):
        result = not self.__eq__(other)
        return result

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
            (self.defined == other.defined) and \
            (abs(self.xdtu - other.xdtu) <= abserror) and \
            (abs(self.ydtu - other.ydtu) <= abserror) and \
            (abs(self.zdtu - other.zdtu) <= abserror) and \
            (abs(self.xdtu_0 - other.xdtu_0) <= abserror) and \
            (abs(self.ydtu_0 - other.ydtu_0) <= abserror) and \
            (abs(self.zdtu_0 - other.zdtu_0) <= abserror)
        return result

    def define_from_fits(self, fitsobj, extnum=0):
        """Define class members from header information in FITS file.

        Parameters
        ----------
        fitsobj: file object
            FITS file whose header contains the DTU information
            needed to initialise the members of this class.
        extnum : int
            Extension number (first extension is 0)

        """

        # read input FITS file
        hdulist = fits.open(fitsobj)
        image_header = hdulist[extnum].header
        hdulist.close()

        # define DTU variables
        list_of_members = self.__dict__.keys()
        list_of_members.remove('defined')
        for member in list_of_members:
            self.__dict__[member] = image_header[member]

        # the attributes have been properly set
        self.defined = True

    def define_from_values(self, xdtu, ydtu, zdtu, xdtu_0, ydtu_0, zdtu_0):
        """Define class members from provided values.

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

        # define DTU variables
        self.xdtu = xdtu
        self.ydtu = ydtu
        self.zdtu = zdtu
        self.xdtu_0 = xdtu_0
        self.ydtu_0 = ydtu_0
        self.zdtu_0 = zdtu_0

        # the attributes have been properly set
        self.defined = True

    def define_from_dictionary(self, inputdict):
        """Define class members from dictionary.

        Parameters
        ----------
        inputdict : dictionary
            Python dictionary defining each member of the class.

        """

        for item in self.__dict__:
            if item != 'defined':
                self.__dict__[item] = inputdict[item]

        # the attributes have been properly set
        self.defined = True

    def outdict(self, precision=3):
        """Return dictionary structure rounded to a given precision."""

        outdict = self.__dict__.copy()
        outdict.pop('defined')
        for item in outdict:
            outdict[item] = round(outdict[item], precision)
        return outdict


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
    list_of_members.remove('defined')

    # compute average of all the members of the class (except 'defined')
    for member in list_of_members:
        result.__dict__[member] = np.mean(
            [tmp_dtu.__dict__[member] for tmp_dtu in list_of_objects]
        )

    result.defined = True

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
    list_of_members.remove('defined')

    # compute maximum difference for each member
    for member in list_of_members:
        tmp_array = np.array(
            [tmp_dtu.__dict__[member] for tmp_dtu in list_of_objects]
        )
        minval = tmp_array.min()
        maxval = tmp_array.max()
        result.__dict__[member] = maxval - minval

    result.defined = True

    return result