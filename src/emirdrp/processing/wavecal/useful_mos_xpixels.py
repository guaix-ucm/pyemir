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

"""Useful X-axis pixels removing +/- npixaround pixels around each OH line"""

import numpy as np
import pkgutil
from six import StringIO

from numina.array.display.ximshow import ximshow
from numina.array.display.pause_debugplot import pause_debugplot

from emirdrp.processing.wavecal.get_islitlet import get_islitlet


def useful_mos_xpixels(reduced_mos_data,
                       base_header,
                       vpix_region,
                       npix_removed_near_ohlines=0,
                       list_valid_wvregions=None,
                       debugplot=0):
    """Useful X-axis pixels removing +/- npixaround pixels around each OH line
    """

    # get wavelength calibration from image header
    naxis1 = base_header['naxis1']
    naxis2 = base_header['naxis2']
    crpix1 = base_header['crpix1']
    crval1 = base_header['crval1']
    cdelt1 = base_header['cdelt1']

    # check vertical region
    nsmin = int(vpix_region[0] + 0.5)
    nsmax = int(vpix_region[1] + 0.5)
    if nsmin > nsmax:
        raise ValueError('vpix_region values in wrong order')
    elif nsmin < 1 or nsmax > naxis2:
        raise ValueError('vpix_region outside valid range')

    # minimum and maximum pixels in the wavelength direction
    islitlet_min = get_islitlet(nsmin)
    ncmin = base_header['jmnslt{:02d}'.format(islitlet_min)]
    ncmax = base_header['jmxslt{:02d}'.format(islitlet_min)]
    islitlet_max = get_islitlet(nsmax)
    if islitlet_max > islitlet_min:
        for islitlet in range(islitlet_min + 1, islitlet_max + 1):
            ncmin_ = base_header['jmnslt{:02d}'.format(islitlet)]
            ncmax_ = base_header['jmnslt{:02d}'.format(islitlet)]
            ncmin = min(ncmin, ncmin_)
            ncmax = max(ncmax, ncmax_)

    # pixels within valid regions
    xisok_wvreg = np.zeros(naxis1, dtype='bool')
    if list_valid_wvregions is None:
        for ipix in range(ncmin, ncmax + 1):
            xisok_wvreg[ipix - 1] = True
    else:
        for wvregion in list_valid_wvregions:
            wvmin = float(wvregion[0])
            wvmax = float(wvregion[1])
            if wvmin > wvmax:
                raise ValueError('wvregion values in wrong order:'
                                 ' {}, {}'.format(wvmin, wvmax))
            minpix = int((wvmin - crval1) / cdelt1 + crpix1 + 0.5)
            maxpix = int((wvmax - crval1) / cdelt1 + crpix1 + 0.5)
            for ipix in range(minpix, maxpix + 1):
                if 1 <= ipix <= naxis1:
                    xisok_wvreg[ipix - 1] = True
        if np.sum(xisok_wvreg) < 1:
            raise ValueError('no valid wavelength ranges provided')

    # pixels affected by OH lines
    xisok_oh = np.ones(naxis1, dtype='bool')
    if int(npix_removed_near_ohlines) > 0:
        dumdata = pkgutil.get_data(
            'emirdrp.instrument.configs',
            'Oliva_etal_2013.dat'
        )
        oh_lines_tmpfile = StringIO(dumdata.decode('utf8'))
        catlines = np.genfromtxt(oh_lines_tmpfile)
        catlines_all_wave = np.concatenate(
            (catlines[:, 1], catlines[:, 0]))
        for waveline in catlines_all_wave:
            expected_pixel = int(
                (waveline - crval1) / cdelt1 + crpix1 + 0.5)
            minpix = expected_pixel - int(npix_removed_near_ohlines)
            maxpix = expected_pixel + int(npix_removed_near_ohlines)
            for ipix in range(minpix, maxpix + 1):
                if 1 <= ipix <= naxis1:
                    xisok_oh[ipix - 1] = False

    # pixels in valid regions not affected by OH lines
    xisok = np.logical_and(xisok_wvreg, xisok_oh)
    naxis1_effective = np.sum(xisok)
    if naxis1_effective < 1:
        raise ValueError('no valid wavelength range available after '
                         'removing OH lines')

    if abs(debugplot) in [21, 22]:
        slitlet2d = reduced_mos_data[(nsmin - 1):nsmax, :].copy()
        ximshow(slitlet2d,
                title='Rectified region',
                first_pixel=(1, nsmin),
                crval1=crval1, cdelt1=cdelt1, debugplot=debugplot)
        ax = ximshow(slitlet2d,
                     title='Rectified region\nafter blocking '
                           'removed wavelength ranges',
                     first_pixel=(1, nsmin),
                     crval1=crval1, cdelt1=cdelt1, show=False)
        for idum in range(1, naxis1 + 1):
            if not xisok_wvreg[idum - 1]:
                ax.plot([idum, idum], [nsmin, nsmax], 'g-')
        pause_debugplot(debugplot, pltshow=True)
        ax = ximshow(slitlet2d,
                     title='Rectified slitlet\nuseful regions after '
                           'removing OH lines',
                     first_pixel=(1, nsmin),
                     crval1=crval1, cdelt1=cdelt1, show=False)
        for idum in range(1, naxis1 + 1):
            if not xisok[idum - 1]:
                ax.plot([idum, idum], [nsmin, nsmax], 'm-')
        pause_debugplot(debugplot, pltshow=True)

    return xisok
