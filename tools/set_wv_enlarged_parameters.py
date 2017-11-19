from __future__ import division
from __future__ import print_function

from emirdrp.core import EMIR_VALID_FILTERS
from emirdrp.core import EMIR_VALID_GRISMS


def set_wv_enlarged_parameters(filter_name, grism_name):
    """Set wavelength calibration parameters for rectified images.

    Parameters
    ----------
    filter_name : str
        Filter name.
    grism_name : str
        Grism name.

    Returns
    -------
    crpix1_enlarged : float
        CRPIX1 value.
    crval1_enlarged : float
        CRVAL1 value.
    cdelt1_enlarged : float
        CDELT1 value.
    naxis1_enlarged: int
        NAXIS1 value.

    """

    # protections
    if filter_name not in EMIR_VALID_FILTERS:
        raise ValueError('Unexpected filter_name:', filter_name)
    if grism_name not in EMIR_VALID_GRISMS:
        raise ValueError('Unexpected grism_name:', grism_name)

    # set parameters
    crpix1_enlarged = 1.0
    if grism_name == "J" and filter_name == "J":
        crval1_enlarged = 11220.0000  # Angstroms
        cdelt1_enlarged = 0.7575      # Angstroms/pixel
        naxis1_enlarged = 3400        # pixels
    elif grism_name == "H" and filter_name == "H":
        crval1_enlarged = None  # 14000.0000  # Angstroms
        cdelt1_enlarged = 1.2000      # Angstroms/pixel
        naxis1_enlarged = 3400        # pixels
    elif grism_name == "K" and filter_name == "Ksp":
        crval1_enlarged = None  # 19000.0000  # Angstroms
        cdelt1_enlarged = 1.7000      # Angstroms/pixel
        naxis1_enlarged = 3400        # pixels
    elif grism_name == "LR" and filter_name == "YJ":
        crval1_enlarged = None        # Angstroms
        cdelt1_enlarged = None        # Angstroms/pixel
        naxis1_enlarged = None        # pixels
    elif grism_name == "LR" and filter_name == "HK":
        crval1_enlarged = None        # Angstroms
        cdelt1_enlarged = None        # Angstroms/pixel
        naxis1_enlarged = None        # pixels
    else:
        print("filter_name..:", filter_name)
        print("grism_name...:", grism_name)
        raise ValueError("invalid filter_name and grism_name combination")

    return crpix1_enlarged, crval1_enlarged, cdelt1_enlarged, naxis1_enlarged
