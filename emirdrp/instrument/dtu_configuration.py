from __future__ import division
from __future__ import print_function

from astropy.io import fits


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
        self.xdtu = image_header['xdtu']
        self.ydtu = image_header['ydtu']
        self.zdtu = image_header['zdtu']
        self.xdtu_0 = image_header['xdtu_0']
        self.ydtu_0 = image_header['ydtu_0']
        self.zdtu_0 = image_header['zdtu_0']

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

    def outdict(self):
        """Return dictionary structure rounded to a given precision."""

        outdict = {}
        if self.defined:
            outdict['xdtu'] = round(self.xdtu, 3)
            outdict['ydtu'] = round(self.ydtu, 3)
            outdict['zdtu'] = round(self.zdtu, 3)
            outdict['xdtu_0'] = round(self.xdtu_0, 3)
            outdict['ydtu_0'] = round(self.ydtu_0, 3)
            outdict['zdtu_0'] = round(self.zdtu_0, 3)

        return outdict
