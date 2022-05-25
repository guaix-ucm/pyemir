#
# Copyright 2015-2021 Universidad Complutense de Madrid
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

import logging
import itertools

import numpy
import scipy.ndimage
from numina.array.utils import expand_region
import numina.array.fwhm as fmod
from numina.array.utils import coor_to_pix_1d
from numina.array.peaks.peakdet import find_peaks_indexes, refine_peaks
from numina.array.utils import slice_create

import emirdrp.instrument.distortions as dist


def find_position(edges, prow, bstart, bend, total=5):
    """Find a EMIR CSU bar position in a edge image.

    Parameters
    ==========
    edges: ndarray,
        a 2d image with 1 where is a border, 0 otherwise
    prow: int,
        reference 'row' of the bars
    bstart: int,
        minimum 'x' position of a bar (0-based)
    bend: int
        maximum 'x' position of a bar (0 based)
    total: int
        number of rows to check near `prow`

    Return
    ======
    list of (x, y) centroids

    """

    nt = total // 2

    # This bar is too near the border
    if prow-nt < 0 or prow + nt >= edges.shape[0]:
        return []

    s2edges = edges[prow-nt:prow+nt+1, bstart:bend]

    structure = scipy.ndimage.generate_binary_structure(2,2) # 8 way conection
    har, num_f = scipy.ndimage.label(s2edges, structure=structure)

    cen_of_mass = scipy.ndimage.center_of_mass(s2edges, labels=har, index=range(1, num_f + 1))

    # center_of_mass returns y, x coordinates

    cen_of_mass_off = [(x + bstart, prow-nt + y) for y,x in cen_of_mass]

    return cen_of_mass_off


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


def simple_prot(x, start):
    """Find the first peak to the right of start"""

    # start must b >= 1

    for i in range(start,len(x)-1):
        a,b,c =  x[i-1], x[i], x[i+1]
        if b - a > 0 and b -c >= 0:
            return i
    else:
        return None


def position_half_h(pslit, cpix, backw=4):
    """Find the position where the value is half of the peak"""

    # Find the first peak to the right of cpix
    next_peak = simple_prot(pslit, cpix)

    if next_peak is None:
        raise ValueError

    dis_peak = next_peak - cpix

    wpos2 = cpix - dis_peak
    wpos1 = wpos2 - backw

    # Compute background in a window of width backw
    # in a position simetrical to the peak
    # around cpix
    left_background = pslit[wpos1:wpos2].min()

    # height of the peak
    height = pslit[next_peak] - left_background

    half_height = left_background + 0.5 * height

    # Position at halg peak, linear interpolation
    vv = pslit[wpos1:next_peak+1] - half_height

    res1, =  numpy.nonzero(numpy.diff(vv > 0))
    i1 = res1[0]

    xint = wpos1 + i1 + (0 - vv[i1]) / (vv[i1+1] - vv[i1])

    return xint, next_peak, wpos1, wpos2, left_background, half_height


def locate_bar_l(icut, epos):
    """Fine position of the left CSU bar"""
    def swap_coor(x):
        return x

    def swap_line(tab):
        return tab

    return _locate_bar_gen(icut, epos,
                           transform1=swap_coor,
                           transform2=swap_line
                           )


def locate_bar_r(icut, epos):
    """Fine position of the right CSU bar"""
    sm = len(icut)

    def swap_coor(x):
        return sm - 1 - x

    def swap_line(tab):
        return tab[::-1]

    return _locate_bar_gen(icut, epos, transform1=swap_coor,
                           transform2=swap_line)


def _locate_bar_gen(icut, epos, transform1, transform2):
    """Generic function for the fine position of the CSU"""

    epos_pix = coor_to_pix_1d(epos)

    # transform ->
    epos_pix_s = transform1(epos_pix)
    icut2 = transform2(icut)
    #

    try:
        res = position_half_h(icut2, epos_pix_s)

        xint_s, next_peak_s, wpos1_s, wpos2_s, background_level, half_height = res
        #

        xint = transform1(xint_s)

        #
        epos_f = xint
        error = 0
    except ValueError:
        error = 2
        epos_f = epos

    return epos_pix, epos_f, error


def char_bar_peak_l(arr_deriv, ypix, bstart, bend, th):
    return _char_bar_peak(arr_deriv, ypix, bstart, bend, th, sign=1)


def char_bar_peak_r(arr_deriv, ypix, bstart, bend, th):
    return _char_bar_peak(arr_deriv, ypix, bstart, bend, th, sign=-1)


def _char_bar_peak(arr_deriv, ypix, bstart, bend, th, sign=1):

    # extract a region to average
    # wy = 3
    # wx = 10
    # Fit the peak with these points
    # wfit = 3

    logger = logging.getLogger(__name__)

    yvpix = numpy.clip(ypix, 0, 2047)

    # Refine at different positions along the slit
    newrefine = []
    wx = 5
    wy = 1
    step = 2 * wy + 1
    # visibility
    #
    intv1 = [0, 2047]
    intv2 = [yvpix - 18, yvpix + 18]
    logger.debug('overlaping interval %s', intv2)

    bar_overlap = overlap(intv1, intv2)
    logger.debug('bar overlaping %.1f', bar_overlap)
    offs = []
    if bar_overlap < 10:
        maxval = 18
    else:
        maxval = 12
    offs2 = range(-step, -maxval, -step)
    offs.extend(reversed(offs2))
    offs.extend(range(0, maxval, step))

    cut = sign * arr_deriv[yvpix, bstart:bend]

    idxs = find_peaks_indexes(cut, window_width=3, threshold=th, fpeak=1)
    logger.debug('found %d peaks over threshold %f', len(idxs), th)

    if len(idxs) == 0:
        logger.debug('no peaks, exit')
        return 0, 0, 0, 0, 0, 1

    # Characterize: use the peak that has the greatest value in the derivative?
    pix_m = cut[idxs].argmax()
    centerx = bstart + idxs[pix_m]
    logger.debug('select the peak with maximum derivative')

    centery = yvpix
    logger.debug('centery is %7.2f at position %7.2f', centery+1,  centerx+1)
    # Refine at the computed center
    xl, fwhm_x, st = refine_bar_centroid(arr_deriv, centerx, centery, wx, wy, th, sign)
    logger.debug('measured values %7.2f (FWHM %7.2f)', xl, fwhm_x)

    if st != 0:
        logger.debug('faillure refining bar centroid, go to next bar')
        # Exiting now, can't refine the centroid
        return centery, centery, xl, xl, fwhm_x, st

    # This is basically to build a list of centers that dont overlap
    for off in offs:
        if 0 <= centery + off <= 2047:
            logger.debug('looping, off %d, measuring at %7.2f', off, centery + off + 1)
            res = refine_bar_centroid(arr_deriv, centerx, centery + off, wx, wy, th, sign)
            logger.debug('looping, measured values %7.2f (FWHM %7.2f)', res[0], res[1])
            newrefine.append(res)
        else:
            logger.debug('looping, off %d, skipping position %7.2f', off, centery + off + 1)

    # this goes in FITS pix coordinates, adding 1
    # filter values with status != 0
    valid_mt = [r[2] == 0 for r in newrefine]
    xcoords_mt = [r[0] + 1 for r in newrefine]
    ycoords_mt = [centery + off + 1 for off in offs]

    ycoords_m = list(itertools.compress(ycoords_mt, valid_mt))
    xcoords_m = list(itertools.compress(xcoords_mt, valid_mt))

    if len(xcoords_m) == 0:
        logger.debug('no valid values to refine')
        return centery, centery, xl, xl, fwhm_x, 3

    logger.debug('transform values from real to virtual')
    xcoords_t, ycoords_t = dist.pvex(xcoords_m, ycoords_m)
    logger.debug('real xcoords are: %s:', xcoords_m)
    logger.debug('real ycoords are: %s:', ycoords_m)
    logger.debug('virtual xcoords are: %s:', xcoords_t)
    logger.debug('virtual ycoords are: %s:', ycoords_t)
    avg_xl_virt = numpy.mean(xcoords_t)
    logger.debug('reference real xcoord is: %s:', xl)
    logger.debug('average virtual xcoord is: %s:', avg_xl_virt)

    centerx_virt, centery_virt = dist.pvex(centerx + 1, centery + 1)

    return centery, centery_virt, xl, avg_xl_virt, fwhm_x, 0


def refine_bar_centroid(arr_deriv, centerx, centery, wx, wy, threshold, sign):
    # Refine values
    logger = logging.getLogger('emir.recipes.bardetect')

    logger.debug('collapsing a %d x %d region', 2 * wx + 1 , 2 * wy + 1)
    #
    slicey = slice_create(centery, wy, start=1, stop=2047)
    slicex = slice_create(centerx, wx, start=1, stop=2047)
    region = arr_deriv[slicey, slicex]

    if region.size == 0:
        logger.debug('region to collapse is empty')
        return 0, 0, 1

    collapsed = sign * region.mean(axis=0)

    # Fine tunning
    idxs_t = find_peaks_indexes(collapsed, window_width=3, threshold=threshold)
    # Use only the peak nearest the original peak
    if len(idxs_t) == 0:
        logger.debug('no peaks after fine-tunning')
        return 0, 0, 2

    dist_t = numpy.abs(idxs_t - wx)
    only_this = dist_t.argmin()
    idxs_p = numpy.atleast_1d(idxs_t[only_this])
    x_t, y_t = refine_peaks(collapsed, idxs_p, window_width=3)

    if len(x_t) == 0:
        logger.debug('no peaks to refine after fitting')
        return 0, 0, 2

    if x_t[0] >= collapsed.shape[0]:
        logger.debug('wrong position %d when refining', x_t[0])
        return 0, 0, 2

    _, fwhm_x = fmod.compute_fwhm_1d_simple(collapsed, x_t[0])

    xl = centerx - wx + x_t[0]
    return xl, fwhm_x, 0


def char_bar_height(arr_deriv_alt, xpos1, xpos2, centery, threshold, wh=35, wfit=3):

    logger = logging.getLogger('emir.recipes.bardetect')
    pcentery = coor_to_pix_1d(centery)
    slicey = slice_create(pcentery, wh, start=1, stop=2047)

    ref_pcentery = pcentery - slicey.start
    mm = arr_deriv_alt[slicey, xpos1:xpos2 + 1].mean(axis=-1)

    idxs_t = find_peaks_indexes(mm, window_width=3,threshold=threshold)
    idxs_u = find_peaks_indexes(-mm, window_width=3,threshold=threshold)
    # Peaks on the right

    status = 0
    npeaks_u = len(idxs_u)
    if npeaks_u == 0:
        # This is a problem, no peak on the right
        b2 = 0
        status = 4
        logger.debug('no bottom border found')
    else:
        # Filter over reference
        g_idxs_u = idxs_u[idxs_u >= ref_pcentery]
        if len(g_idxs_u) == 0:
            logger.debug('no peak over center')
            b2 = 0
            status = 4
        else:
            x_u, y_u = refine_peaks(-mm, g_idxs_u, window_width=wfit)
            # Select the peak with max derivative
            if len(x_u) == 0 or len(y_u) == 0:
                logger.warning("no 1st peak found after refine")
                b2 = 0
                status = 4
            else:
                idmax = y_u.argmax()
                b2 = x_u[idmax]
                b2val = y_u[idmax]
                logger.debug('main border in %f', slicey.start + b2)

    # peaks on the left
    npeaks_t = len(idxs_t)
    if npeaks_t == 0:
        # This is a problem, no peak on the left
        b1 = 0
        logger.debug('no top border found')
        status = 40 + status
    else:
        g_idxs_t = idxs_t[idxs_t <= ref_pcentery]
        if len(g_idxs_t) == 0:
            logger.debug('no peak under center')
            b1 = 0
            status = 40 + status
        else:
            x_t, y_t = refine_peaks(mm, g_idxs_t, window_width=wfit)
            # Select the peak with max derivative
            if len(x_t) == 0 or len(y_t) == 0:
                logger.warning("no 2nd peak found after refine")
                b1 = 0
                status = 40 + status
            else:
                idmax = y_t.argmax()
                b1 = x_t[idmax]
                b1val = y_t[idmax]
                logger.debug('second border in %f', slicey.start + b1)

    return slicey.start + b1, slicey.start + b2, status


def overlap(intv1, intv2):
    """Overlaping of two intervals"""
    return max(0, min(intv1[1], intv2[1]) - max(intv1[0], intv2[0]))