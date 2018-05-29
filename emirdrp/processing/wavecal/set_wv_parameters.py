#
# Copyright 2018 Universidad Complutense de Madrid
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

"""Wavelength calibration parameters for each grism + filter configuration"""

from __future__ import division, print_function

import numpy as np

from emirdrp.core import EMIR_NAXIS1_ENLARGED
from emirdrp.core import EMIR_VALID_FILTERS
from emirdrp.core import EMIR_VALID_GRISMS


def set_wv_parameters(filter_name, grism_name):
    """Set wavelength calibration parameters for rectified images.

    Parameters
    ----------
    filter_name : str
        Filter name.
    grism_name : str
        Grism name.

    Returns
    -------
    wv_parameters : dictionary
        Python dictionary containing relevant wavelength calibration
        parameters:
        - crpix1_enlarged : float
        - crval1_enlarged : float
        - cdelt1_enlarged : float
        - naxis1_enlarged: int
        - islitlet_min : int
          Minimium slitlet number.
        - islitlet_max : int
          Maximium slitlet number.
        - nbrightlines : python list
          List of integers containing the number of brightlines to be
          used in the initial wavelength calibration.
        - poly_crval1_linear : numpy.polynomial.Polynomial instance
          Polynomial providing the value of CRVAL1_linear as a function
          of csu_bar_slit_center.
        - poly_cdelt1_linear : numpy.polynomial.Polynomial instance
          Polynomial providing the value of CDELT1_linear as a function
          of csu_bar_slit_center.

    """

    # protections
    if filter_name not in EMIR_VALID_FILTERS:
        raise ValueError('Unexpected filter_name:', filter_name)
    if grism_name not in EMIR_VALID_GRISMS:
        raise ValueError('Unexpected grism_name:', grism_name)

    # intialize output
    wv_parameters = {}
    # set parameters
    wv_parameters['crpix1_enlarged'] = 1.0
    if grism_name == "J" and filter_name == "J":
        wv_parameters['islitlet_min'] = 2
        wv_parameters['islitlet_max'] = 54
        wv_parameters['nbrightlines'] = [18]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            1.25122231e+04,
            -4.83412417e+00,
            5.31657015e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            7.73411692e-01,
            -3.28055653e-05,
            -1.20813896e-08
        ])
        wv_parameters['crval1_enlarged'] = 11200.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 0.77        # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = EMIR_NAXIS1_ENLARGED  # pixels
        wv_parameters['wvmin_expected'] = None
        wv_parameters['wvmax_expected'] = None
        wv_parameters['wvmin_useful'] = None
        wv_parameters['wvmax_useful'] = None
    elif grism_name == "H" and filter_name == "H":
        wv_parameters['islitlet_min'] = 2
        wv_parameters['islitlet_max'] = 54
        wv_parameters['nbrightlines'] = [12]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            1.65536274e+04,
            -7.63517173e+00,
            7.74790265e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            1.21327515e+00,
            1.42140078e-05,
            -1.27489119e-07
        ])
        wv_parameters['crval1_enlarged'] = 14500.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 1.2200      # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = EMIR_NAXIS1_ENLARGED  # pixels
        wv_parameters['wvmin_expected'] = None
        wv_parameters['wvmax_expected'] = None
        wv_parameters['wvmin_useful'] = None
        wv_parameters['wvmax_useful'] = None
    elif grism_name == "K" and filter_name == "Ksp":
        wv_parameters['islitlet_min'] = 2
        wv_parameters['islitlet_max'] = 54
        wv_parameters['nbrightlines'] = [12]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            2.21044741e+04,
            -1.08737529e+01,
            9.05081653e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            1.72696857e+00,
            2.35009351e-05,
            -1.02164228e-07
        ])
        wv_parameters['crval1_enlarged'] = 19100.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 1.7300      # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = EMIR_NAXIS1_ENLARGED  # pixels
        wv_parameters['wvmin_expected'] = None
        wv_parameters['wvmax_expected'] = None
        wv_parameters['wvmin_useful'] = None
        wv_parameters['wvmax_useful'] = None
    elif grism_name == "LR" and filter_name == "YJ":
        wv_parameters['islitlet_min'] = 4
        wv_parameters['islitlet_max'] = 55
        wv_parameters['nbrightlines'] = [20]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            1.04272465e+04,
            -2.33176855e+01,
            6.55101267e-03
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            3.49037727e+00,
            1.26008332e-03,
            -4.66149678e-06
        ])
        wv_parameters['crval1_enlarged'] = 8900.0000   # Angstroms
        wv_parameters['cdelt1_enlarged'] = 3.5600      # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = 1270        # pixels
        wv_parameters['wvmin_expected'] = 4000
        wv_parameters['wvmax_expected'] = 15000
        wv_parameters['wvmin_useful'] = 8900
        wv_parameters['wvmax_useful'] = 13400
    elif grism_name == "LR" and filter_name == "HK":
        wv_parameters['islitlet_min'] = 4
        wv_parameters['islitlet_max'] = 55
        wv_parameters['nbrightlines'] = [25]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            2.00704978e+04,
            -4.07702886e+01,
            -5.95247468e-03
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            6.54247758e+00,
            2.09061196e-03,
            -2.48206609e-06
        ])
        wv_parameters['crval1_enlarged'] = 14500.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 6.83        # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = 1435        # pixels
        wv_parameters['wvmin_expected'] = 8000
        wv_parameters['wvmax_expected'] = 30000
        wv_parameters['wvmin_useful'] = 14500
        wv_parameters['wvmax_useful'] = 24300
    else:
        print("filter_name..:", filter_name)
        print("grism_name...:", grism_name)
        raise ValueError("invalid filter_name and grism_name combination")

    return wv_parameters
