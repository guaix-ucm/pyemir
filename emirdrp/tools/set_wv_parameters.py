from __future__ import division
from __future__ import print_function

import numpy as np

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
            1.25137094e+04,
            -4.81553731e+00,
            4.70039758e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            7.74133267e-01,
            -4.72423718e-05,
            2.79842624e-08
        ])
        wv_parameters['crval1_enlarged'] = 11200.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 0.77        # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = 3400        # pixels
    elif grism_name == "H" and filter_name == "H":
        wv_parameters['islitlet_min'] = 2
        wv_parameters['islitlet_max'] = 54
        wv_parameters['nbrightlines'] = [15]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            1.65561291e+04,
            -7.62779414e+00,
            7.52228172e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            1.21372302e+00,
            3.77209658e-06,
            -1.07467162e-07
        ])
        wv_parameters['crval1_enlarged']= 14500.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 1.2200      # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = 3400        # pixels
    elif grism_name == "K" and filter_name == "Ksp":
        wv_parameters['islitlet_min'] = 2
        wv_parameters['islitlet_max'] = 54
        wv_parameters['nbrightlines'] = [15]
        wv_parameters['poly_crval1_linear'] = np.polynomial.Polynomial([
            2.21095313e+04,
            -1.08900414e+01,
            9.66474839e-04
        ])
        wv_parameters['poly_cdelt1_linear'] = np.polynomial.Polynomial([
            1.72596244e+00,
            2.85573046e-05,
            -1.30027272e-07
        ])
        wv_parameters['crval1_enlarged'] = 19100.0000  # Angstroms
        wv_parameters['cdelt1_enlarged'] = 1.7300      # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = 3400        # pixels
    elif grism_name == "LR" and filter_name == "YJ":
        wv_parameters['islitlet_min'] = 4
        wv_parameters['islitlet_max'] = 55
        wv_parameters['nbrightlines'] = None
        wv_parameters['poly_crval1_linear'] = None
        wv_parameters['poly_cdelt1_linear'] = None
        wv_parameters['crval1_enlarged'] = None        # Angstroms
        wv_parameters['cdelt1_enlarged'] = None        # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = None        # pixels
    elif grism_name == "LR" and filter_name == "HK":
        wv_parameters['islitlet_min'] = 4
        wv_parameters['islitlet_max'] = 55
        wv_parameters['nbrightlines'] = None
        wv_parameters['poly_crval1_linear'] = None
        wv_parameters['poly_cdelt1_linear'] = None
        wv_parameters['crval1_enlarged'] = None        # Angstroms
        wv_parameters['cdelt1_enlarged'] = None        # Angstroms/pixel
        wv_parameters['naxis1_enlarged'] = None        # pixels
    else:
        print("filter_name..:", filter_name)
        print("grism_name...:", grism_name)
        raise ValueError("invalid filter_name and grism_name combination")

    return wv_parameters
