#
# Copyright 2015 Universidad Complutense de Madrid
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

"""Bar detection procedures for EMIR"""

import numpy
from scipy.ndimage.filters import median_filter
from skimage.feature import canny

from numina.array.utils import expand_region
import numina.array.fwhm as fmod
from numina.array.utils import wc_to_pix_1d

from .common import normalize_raw


def load_position_table():
    fname = "final.csv"
    table = numpy.loadtxt(fname, delimiter=',')
    # Coordinates are in FITS, so we
    # substract 1
    return table - 1


def find_position(edges, yref, bstart, bend, total=5, maxdist=1.5):
    """Find a EMIR CSU bar position in a edge image.

    Parameters
    ==========
    edges; ndarray,
        a 2d image with 1 where is a border, 0 otherwise
    yref: float,
        reference 'y' coordinate of this bar
    bstart: int,
        minimum 'x' position of a bar (0-based)
    bend: int
        maximum 'x' position of a bar (0 based)
    total: int
        number of rows to check near `yref`
    maxdist: float
        maximum distance between peaks in different rows

    Return
    ======
    tuple with row as given by yref, left border and rigth border of the bar,
    None if the bar is not found

    """
    prow = wc_to_pix_1d(yref)

    nt = total // 2

    cents = []
    # do "total" cuts and find peaks
    for h in range(-nt, nt+1):
        sedges = edges[prow+h, bstart:bend]
        cuts, = numpy.nonzero(sedges==1)
        tcuts = cuts + bstart
        # if there are exactly 2 peaks
        # accumulate these pair of borders
        if len(tcuts) > 2:
            continue
        if len(tcuts) < 2:
            continue

        cents.append(tcuts)

    ncents = numpy.array(cents)

    # skip this array, is empty
    # we can't find a bar here
    if ncents.ndim != 2:
        return None

    # find the mean of positions of peaks
    # if the distance to the reference
    # is less than maxdist
    m = numpy.abs(ncents - cents[nt])
    fc = m[:,0] < maxdist
    fd = m[:,1] < maxdist

    c1 = ncents[fc,0].mean()
    c2 = ncents[fd,1].mean()

    thisres = (prow, c1, c2)

    return thisres


def calc_fwhm(img, region, fexpand=3, axis=0):
    """Compute the FWHM in the direction given by axis"""

    # We compute know the FWHM of the slit
    # Given the computed position of the slit
    # Expand 'fexpand' pixels around
    # and cut an slice in the median filtered image

    xpregion = expand_region(region, fexpand, fexpand)
    cslit = img[xpregion]

    # Collapse it
    pslit = cslit.mean(axis=axis)

    # Estimate the background as a flat line
    # starting in pslit[0] and ending in pslit[-1]
    x2 = len(pslit)
    y1, y2 = pslit[0], pslit[-1]
    mslope = (y2-y1) / x2
    # background estimation
    backstim = mslope*numpy.arange(x2) + y1

    # We subtract background
    qslit = pslit-backstim
    # and find the pixel of the maximum
    pidx = numpy.argmax(qslit)
    peak, fwhm = fmod.compute_fwhm_1d_simple(qslit, pidx)
    return fwhm


def recipe_function(arr, slitstab):

    # Median filter
    mfilter_size = 5

    arr_median = median_filter(arr, size=mfilter_size)

    # Image is mapped between 0 and 1
    # for the full range [0: 2**16]
    arr_grey = normalize_raw(arr_median)

    # Find borders
    canny_sigma = 3.0
    # These threshols corespond roughly with
    # value x (2**16 - 1)
    high_threshold = 0.04
    low_threshold = 0.01
    edges = canny(arr_grey, sigma=canny_sigma,
                  high_threshold=high_threshold,
                  low_threshold=low_threshold)

    # Number or rows used

    total = 5
    maxdist = 1.0
    nt = total // 2

    bstart = 100
    bend = 1900

    fexpand = 3

    result = []

    # Based om the 'edges image'
    # and the table of approx positions of the slits
    for slitid, coords in enumerate(slitstab):

        # Find the position of each bar
        bpos = find_position(edges, coords[1], bstart, bend, total, maxdist)

        # If no bar is found, append and empty token
        if bpos is None:
            thisres = (slitid, )
        else:
            prow, c1, c2 = bpos

            # Compute FWHM of the collapsed profile

            region = (slice(prow-nt, prow+nt+1), slice(c1, c2+1))
            fwhm = calc_fwhm(arr_grey, region, fexpand)

            thisres = (slitid, prow+1, c1+1, c2+1, fwhm)

        result.append(thisres)


    return result
