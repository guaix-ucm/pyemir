import numpy
import astropy.units as u
import astropy.wcs
from astropy import wcs as wcs

import emirdrp.instrument.constants as cons


def create_wcs(fppa=270 * u.deg):

    w = astropy.wcs.WCS(naxis=2)
    ipa = cons.EMIR_REF_IPA
    angle = fppa - ipa
    # scale = numpy.cos(numpy.deg2rad(ipa - fppa))
    w.wcs.ctype = ["RA---ZPN", "DEC--ZPN"]
    w.wcs.crval = [86.3712733276122, 28.9182815703162]
    w.wcs.crpix = [1022.68, 1016.70]
    # scale = cons.EMIR_PIXSCALE.to(u.deg / u.pixel).value
    scale1 = 0.194279723 / 3600
    scale2 = 0.194264445 / 3600
    # w.wcs.cdelt = [scale, scale]
    pc11 = scale1 * numpy.cos(angle)
    pc22 = scale2 * numpy.cos(angle)
    pc12 = -scale1 * numpy.sin(angle)
    pc21 = scale2 * numpy.sin(angle)
    w.wcs.cd = [[pc11, pc12], [pc21, pc22]]
    w.wcs.set_pv(
        [(2, 1, 0.9998), (2, 2, 0.0), (2, 3, 13824.76), (2, 4, 0.0), (2, 5, 3491446467)]
    )
    w.wcs.radesys = "FK5"
    return w


def create_wcs_new():
    """From EMIR image as 2020-07-05"""
    w = astropy.wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---ZPN", "DEC--ZPN"]
    w.wcs.crval = [297.591888896875, 30.7466139938489]
    w.wcs.crpix = [1022.46378968092, 1016.55073697444]
    pc11 = -5.39644538444764e-05
    pc22 = -5.39644538444764e-05
    pc12 = -5.19905444897492e-08
    pc21 = 5.19905444897492e-08
    w.wcs.cd = [[pc11, pc12], [pc21, pc22]]
    w.wcs.set_pv(
        [
            (2, 1, 1),
            (2, 2, 0.0),
            (2, 3, 14584.8),
            (2, 4, 6.93534705104727e-310),
            (2, 5, 1755295542.4),
        ]
    )
    w.wcs.radesys = "FK5"
    return w


def create_wcs_alt(fppa=270 * u.deg):
    """Candidate alternate WCS to pass to virtual pixels, not working"""

    w = astropy.wcs.WCS(naxis=2)
    ipa = cons.EMIR_REF_IPA
    angle = fppa - ipa
    # scale = numpy.cos(numpy.deg2rad(ipa - fppa))
    w.wcs.ctype = ["", ""]
    w.wcs.crval = [86.3712733276122, 28.9182815703162]
    w.wcs.crpix = [1022.68, 1016.70]
    # scale = cons.EMIR_PIXSCALE.to(u.deg / u.pixel).value
    scale1 = 0.194279723 / 3600
    scale2 = 0.194264445 / 3600
    # w.wcs.cdelt = [scale, scale]
    pc11 = scale1 * numpy.cos(angle)
    pc22 = scale2 * numpy.cos(angle)
    pc12 = -scale1 * numpy.sin(angle)
    pc21 = scale2 * numpy.sin(angle)
    w.wcs.cd = [[pc11, pc12], [pc21, pc22]]
    w.wcs.radesys = "FK5"
    return w


def create_wcs_4829(crpix=None):

    w = wcs.WCS(naxis=2)

    crpix = [50.0, 50.0] if crpix is None else crpix

    w.wcs.crpix = crpix
    w.wcs.crval = [0.0, 0.0]
    w.wcs.cdelt = [1.0, 1.0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cd = numpy.array(
        [[5.41615466e-05, -2.55232165e-07], [2.55232165e-07, 5.41615466e-05]]
    )
    return w
