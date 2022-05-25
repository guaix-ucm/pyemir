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

import numpy as np
from astropy.io import fits


def save_ndarray_to_fits(array=None, file_name=None,
                         main_header=None,
                         cast_to_float=True,
                         crpix1=None, crval1=None, cdelt1=None,
                         overwrite=True):
    """Save numpy array(s) into a FITS file with the provided filename.

    Parameters
    ----------
    array : numpy array or list of numpy arrays
        Array(s) to be exported as the FITS file. If the input is
        a list, a multi-extension FITS file is generated assuming
        that the list contains a list of arrays.
    file_name : string
        File name for the output FITS file.
    main_header : astropy FITS header
        Header to be introduced in the primary HDU.
    cast_to_float : bool or list of booleans
        If True, the array(s) data are save as float. If a list of
        arrays has been provided, this parameter must be either a list
        (with the same length) of booleans or None.
    crpix1 : float, list of floats or None
        If not None, this value is used for the keyword CRPIX1.
        If a list of arrays has been provided, this paramater must be
        either a list (with the same length) of floats or None.
    crval1 : float, list of floats or None
        If not None, this value is used for the keyword CRVAL1.
        If a list of arrays has been provided, this paramater must be
        either a list (with the same length) of floats or None.
    cdelt1 : float, list of floats or None
        If not None, this value is used for the keyword CDELT1.
        If a list of arrays has been provided, this paramater must be
        either a list (with the same length) of floats or None.
    overwrite : bool
        If True, the file is overwritten (in the case it already
        exists).

    """

    # protections
    if file_name is None:
        raise ValueError("File_name is not defined in save_ndarray_to_fits")

    if type(array) is list:
        list_of_arrays = array
        narrays = len(list_of_arrays)
        # cast_to_float must be a list of bools
        if type(cast_to_float) is not list:
            raise ValueError("Expected list of cast_to_float not found!")
        else:
            if len(cast_to_float) != narrays:
                raise ValueError("Unexpected length of cast_to_float")
        list_cast_to_float = cast_to_float
        # check that the additional associated lists have been provided
        # and that they have the expected length (or they are None)
        for ldum, cdum in zip([crpix1, crval1, cdelt1],
                              ['crpix1', 'crval1', 'cdelt1']):
            if ldum is not None:
                if type(ldum) is not list:
                    raise ValueError("Expected list of " + cdum +
                                     " not found!")
                else:
                    if len(ldum) != narrays:
                        raise ValueError("Unexpected length of " + cdum)
        if crpix1 is None:
            list_crpix1 = [None] * narrays
        else:
            list_crpix1 = crpix1
        if crval1 is None:
            list_crval1 = [None] * narrays
        else:
            list_crval1 = crval1
        if cdelt1 is None:
            list_cdelt1 = [None] * narrays
        else:
            list_cdelt1 = cdelt1
    else:
        list_of_arrays = [array]
        list_cast_to_float = [cast_to_float]
        list_crpix1 = [crpix1]
        list_crval1 = [crval1]
        list_cdelt1 = [cdelt1]

    hdulist = fits.HDUList()

    for ihdu, tmp_array in enumerate(list_of_arrays):
        if type(tmp_array) is not np.ndarray:
            raise ValueError("Array#" + str(ihdu) + "=" + str(tmp_array) +
                             " must be a numpy.ndarray")
        if ihdu == 0:
            if list_cast_to_float[ihdu]:
                hdu = fits.PrimaryHDU(data=tmp_array.astype(np.float32),
                                      header=main_header)
            else:
                hdu = fits.PrimaryHDU(data=tmp_array, header=main_header)
        else:
            if list_cast_to_float[ihdu]:
                hdu = fits.ImageHDU(data=tmp_array.astype(np.float32))
            else:
                hdu = fits.ImageHDU(data=tmp_array)

        # set additional FITS keywords if requested
        tmp_crpix1 = list_crpix1[ihdu]
        if tmp_crpix1 is not None:
            hdu.header.set('CRPIX1', tmp_crpix1, 'Reference pixel')
        tmp_crval1 = list_crval1[ihdu]
        if tmp_crval1 is not None:
            hdu.header.set('CRVAL1', tmp_crval1,
                           'Reference wavelength corresponding to CRPIX1')
        tmp_cdelt1 = list_cdelt1[ihdu]
        if tmp_cdelt1 is not None:
            hdu.header.set('CDELT1', tmp_cdelt1,
                           'Linear dispersion (angstrom/pixel)')

        # add HDU to HDUList
        hdulist.append(hdu)

    # write output file
    hdulist.writeto(file_name, overwrite=overwrite)
