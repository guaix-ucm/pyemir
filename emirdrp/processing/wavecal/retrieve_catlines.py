#
# Copyright 2019 Universidad Complutense de Madrid
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

import numpy as np
import pkgutil
from six import StringIO


def retrieve_catlines(mode, grism_name):
    """Retrieve arc/OH lines

    Parameters
    ----------
    mode : int
        Integer, indicating the type of lines:
        1 or 2 : arc lines
        11 or 12 : OH lines
    grism_name : string
        Grism name.

    Returns
    -------
    catlines_all_wave : numpy array
        Array with wavelengths
    catlines_all_flux : numpy array
        Array with fluxes.

    """

    # protections
    if mode not in [1, 2, 11, 12]:
        raise ValueError('Invalid mode={}'.format(mode))

    # read tabulated lines
    if mode in [1, 2]:  # ARC lines
        if grism_name == 'LR':
            catlines_file = 'lines_argon_neon_xenon_empirical_LR.dat'
        else:
            catlines_file = 'lines_argon_neon_xenon_empirical.dat'
        dumdata = pkgutil.get_data('emirdrp.instrument.configs', catlines_file)
        arc_lines_tmpfile = StringIO(dumdata.decode('utf8'))
        catlines = np.genfromtxt(arc_lines_tmpfile)
        # define wavelength and flux as separate arrays
        catlines_all_wave = catlines[:, 0]
        catlines_all_flux = catlines[:, 1]
    elif mode in [11, 12]:  # OH lines
        dumdata = pkgutil.get_data(
            'emirdrp.instrument.configs',
            'Oliva_etal_2013.dat'
        )
        oh_lines_tmpfile = StringIO(dumdata.decode('utf8'))
        catlines = np.genfromtxt(oh_lines_tmpfile)
        # define wavelength and flux as separate arrays
        catlines_all_wave = np.concatenate((catlines[:, 1], catlines[:, 0]))
        catlines_all_flux = np.concatenate((catlines[:, 2], catlines[:, 2]))
    else:
        raise ValueError('Unexpected mode={}'.format(mode))

    return catlines_all_wave, catlines_all_flux
