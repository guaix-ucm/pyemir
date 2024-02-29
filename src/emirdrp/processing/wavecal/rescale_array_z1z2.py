#
# Copyright 2018-2024 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Auxiliary functions to rescale numpy arrays"""

import numpy as np


def rescale_array_to_z1z2(array, z1z2=(-1.0, 1.0)):
    """Rescale the values in a numpy array to the [z1,z2] interval.

    The transformation is carried out following the relation
    array_rs = b_flux * array - c_flux
    as explained in Appendix B1 of Cardiel (2009, MNRAS, 396, 680)

    Parameters
    ----------
    array : numpy array
        Numpy array to be rescaled.
    z1z2 : tuple, floats
        Minimum and maximum values in the returned array.

    Returns
    -------
    array_rs : numpy array
        Array with rescaled values.
    coef_rs : tuple, floats
        Coefficients b_flux and c_flux employed in the rescaling
        operation.

    """

    if type(array) is not np.ndarray:
        raise ValueError("array=" + str(array) + " must be a numpy.ndarray")

    array_min = array.min()
    array_max = array.max()

    z1, z2 = z1z2
    delta = array_max - array_min
    b_flux = (z2 - z1) / delta
    c_flux = (z2 * array_min - z1 * array_max) / delta

    array_rs = b_flux * array - c_flux

    return array_rs, (b_flux, c_flux)


def rescale_array_from_z1z2(array_rs, coef_rs=None):
    """Restore the values in a numpy array rescaled to the [z1,z2] interval.

    The transformation is carried out following the relation
    array = (array_rs + c_flux)/b_flux
    as explained in Appendix B1 of Cardiel (2009, MNRAS, 396, 680)

    Parameters
    ----------
    array_rs : numpy array
        Numpy array previously rescaled to the [z1,z2] interval with
        the function rescale_array_to_z1z2().
    coef_rs : tuple, floats
        Coefficients b_flux and c_flux previously employed in the
        rescaling operation. This tuple is one of the parameters
        returned by function_rescale_array_to_z1z2().

    Returns
    -------
    array : numpy array
        Array with restored values.

    """

    if type(array_rs) is not np.ndarray:
        raise ValueError(
            "array_rs=" + str(array_rs) + "must be a numpy.ndarray")

    b_flux, c_flux = coef_rs

    array = (array_rs + c_flux) / b_flux

    return array
