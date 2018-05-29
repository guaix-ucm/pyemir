#
# Copyright 2013-2018 Universidad Complutense de Madrid
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

"""AIV Recipes for EMIR"""

from __future__ import division

import logging
import math

import numpy
import scipy.interpolate as itpl
import scipy.optimize as opz
from scipy import ndimage

from astropy.modeling import models, fitting
import astropy.wcs
import photutils
from photutils import CircularAperture

from numina.array.recenter import centering_centroid
from numina.array.utils import image_box
from numina.array.fwhm import compute_fwhm_2d_spline
from numina.constants import FWHM_G
from numina.array.fwhm import compute_fwhm_2d_simple
from numina.array.utils import expand_region

from emirdrp.core import EMIR_PIXSCALE
from .procedures import compute_fwhm_enclosed_direct
from .procedures import compute_fwhm_enclosed_grow
from .procedures import moments
from .procedures import AnnulusBackgroundEstimator
from .procedures import image_box2d

_logger = logging.getLogger(__name__)


# returns y,x
def compute_fwhm(img, center):
    X = numpy.arange(img.shape[0])
    Y = numpy.arange(img.shape[1])

    bb = itpl.RectBivariateSpline(X, Y, img)
    # We assume that the peak is in the center...
    peak = bb(*center)[0, 0]

    def f1(x):
        return bb(x, center[1]) - 0.5 * peak

    def f2(y):
        return bb(center[0], y) - 0.5 * peak

    def compute_fwhm_1(U, V, fun, center):

        cp = int(math.floor(center + 0.5))

        # Min on the rigth
        r_idx = V[cp:].argmin()
        u_r = U[cp + r_idx]

        if V[cp + r_idx] > 0.5 * peak:
            # FIXME: we have a problem
            # brentq will raise anyway
            pass

        sol_r = opz.brentq(fun, center, u_r)

        # Min in the left
        rV = V[cp - 1::-1]
        rU = U[cp - 1::-1]
        l_idx = rV.argmin()
        u_l = rU[l_idx]
        if rV[l_idx] > 0.5 * peak:
            # FIXME: we have a problem
            # brentq will raise anyway
            pass

        sol_l = opz.brentq(fun, u_l, center)
        fwhm = sol_r - sol_l
        return fwhm

    U = X
    V = bb.ev(U, [center[1] for _ in U])
    fwhm_x = compute_fwhm_1(U, V, f1, center[0])

    U = Y
    V = bb.ev([center[0] for _ in U], U)

    fwhm_y = compute_fwhm_1(U, V, f2, center[1])

    return center[0], center[1], peak, fwhm_x, fwhm_y


# Background in an annulus, mode is HSM
def compute_fwhm_global(data, center, box):
    sl = image_box(center, data.shape, box)
    raster = data[sl]

    background = raster.min()
    braster = raster - background

    newc = center[0] - sl[0].start, center[1] - sl[1].start
    try:
        res = compute_fwhm(braster, newc)
        return (res[1] + sl[1].start, res[0] + sl[0].start,
                res[2], res[3], res[4])
    except ValueError as error:
        _logger.warning("%s", error)
        return center[1], center[0], -99.0, -99.0, -99.0
    except Exception as error:
        _logger.warning("%s", error)
        return center[1], center[0], -199.0, -199.0, -199.0


# returns x,y
def gauss_model(data, center_r):
    sl = image_box(center_r, data.shape, box=(4, 4))
    raster = data[sl]

    # background
    background = raster.min()

    b_raster = raster - background

    new_c = center_r[0] - sl[0].start, center_r[1] - sl[1].start

    yi, xi = numpy.indices(b_raster.shape)

    g = models.Gaussian2D(amplitude=b_raster.max(), x_mean=new_c[
                          1], y_mean=new_c[0], x_stddev=1.0, y_stddev=1.0)
    f1 = fitting.LevMarLSQFitter()  # @UndefinedVariable
    t = f1(g, xi, yi, b_raster)

    mm = t.x_mean.value + sl[1].start, t.y_mean.value + \
        sl[0].start, t.amplitude.value, t.x_stddev.value, t.y_stddev.value
    return mm

 
def recenter_char(data, centers_i, recenter_maxdist, recenter_nloop, recenter_half_box, do_recenter):
    
    # recentered values
    centers_r = numpy.empty_like(centers_i)
    # Ignore certain pinholes
    compute_mask = numpy.ones((centers_i.shape[0],), dtype='bool')
    status_array = numpy.ones((centers_i.shape[0],), dtype='int')
    
    for idx, (xi, yi) in enumerate(centers_i):
        # A failsafe
        _logger.info('for pinhole %i', idx)
        _logger.info('center is x=%7.2f y=%7.2f', xi, yi)
        if (xi > data.shape[1] - 5 or xi < 5 or
                yi > data.shape[0] - 5 or yi < 5):
            _logger.info('pinhole too near to the border')
            compute_mask[idx] = False
            centers_r[idx] = xi, yi
            status_array[idx] = 0
        else:
            if do_recenter and (recenter_maxdist > 0.0):
                kk = centering_centroid(
                    data, xi, yi,
                    box=recenter_half_box,
                    maxdist=recenter_maxdist, nloop=recenter_nloop
                )
                xc, yc, _back, status, msg = centering_centroid(
                    data, xi, yi,
                    box=recenter_half_box,
                    maxdist=recenter_maxdist, nloop=recenter_nloop
                )
                _logger.info('new center is x=%7.2f y=%7.2f', xc, yc)
                # Log in X,Y format
                _logger.debug('recenter message: %s', msg)
                centers_r[idx] = xc, yc
                status_array[idx] = status
            else:
                centers_r[idx] = xi, yi
                status_array[idx] = 0

    return centers_r, compute_mask, status_array


def pinhole_char(data, ncenters, box=4, recenter_pinhole=True, maxdist=10.0):

    ibox = (box, box)

    # convert FITS x, y coordinates (pixel center 1)
    # to python/scipy/astropy (pixel center 0 and x,y -> y,x)

    centers_i = ncenters[:, 0:2] - 1
    centers_r, cmask, starr = recenter_char(data, centers_i,
                                            recenter_maxdist=maxdist,
                                            recenter_nloop=10,
                                            recenter_half_box=ibox,
                                            do_recenter=recenter_pinhole
                                            )

    # Result 0
    nresults = 11
    mm0 = numpy.empty((centers_r.shape[0], nresults))
    mm0[:, 0:2] = centers_r
    mm0[:, 2] = starr
    mm0[:, 3:] = -99

    # compute the FWHM without fitting
    for idx, (xc, yc) in enumerate(centers_r):
        _logger.info('For pinhole %i', idx)
        if cmask[idx]:
            _logger.info('compute model-free FWHM')
            try:
                res1 = compute_fwhm_global(data, (yc, xc), box=ibox)
                fmt1 = 'x=%7.2f y=%7.2f peak=%6.3f fwhm_x=%6.3f fwhm_y=%6.3f'
                _logger.info(fmt1, *res1)
                mm0[idx, 3:6] = res1[2:]
            except Exception as error:
                _logger.exception("unable to obtain FWHM, %s", error)

            _logger.info('compute Gaussian 2Dfitting')
            try:
                res2 = gauss_model(data, (yc, xc))
                fmt2 = 'x=%7.2f y=%7.2f peak=%6.3f stdev_x=%6.3f stdev_y=%6.3f'
                _logger.info(fmt2, *res2)
                mm0[idx, 6:9] = res2[2:]
            except Exception as error:
                _logger.exception("unable to obtain FWHM, %s", error)
        else:
            _logger.info('skipping')

    # Photometry in coordinates
    # x=centers_r[:,0]
    # y=centers_r[:,1]
    # with radius 2.0 pixels and 4.0 pixels

    apertures = [2.0, 4.0]
    for idx, rad in enumerate(apertures):
        _logger.info('compute photometry with aperture radii %s', rad)
        apertures = CircularAperture(centers_r[:,:2], r=rad)
        faper = photutils.aperture_photometry(data, apertures)
        mm0[:, 9+idx] = faper['aperture_sum']
        _logger.info('done')

    # Convert coordinates to FITS
    mm0[:, 0:2] += 1
    return mm0


def pinhole_char2(
    data, ncenters,
    recenter_pinhole=True,
    recenter_half_box=5,
    recenter_nloop=10,
    recenter_maxdist=10.0,
    back_buff=3,
    back_width=5,
    phot_niter=10,
    phot_rad=8
):

    sigma0 = 1.0
    rad = 3 * sigma0 * FWHM_G
    box = recenter_half_box
    recenter_half_box = (box, box)

    # convert FITS x, y coordinates (pixel center 1)
    # to python/scipy/astropy (pixel center 0 and x,y -> y,x)

    centers_i = ncenters[:, 0:2] - 1

    centers_r, cmask, starr = recenter_char(data, centers_i, recenter_maxdist,
                                            recenter_nloop, recenter_half_box,
                                            do_recenter=recenter_pinhole
                                            )

    # Number of results
    nresults = 35
    mm0 = numpy.empty((centers_r.shape[0], nresults))
    mm0[:, 0:2] = centers_i
    mm0[:, 2:4] = centers_r
    mm0[:, 4] = starr
    mm0[:, 5:] = -99

    # Fitter
    fitter = fitting.LevMarLSQFitter()  # @UndefinedVariable
    rplot = phot_rad

    for idx, (x0, y0) in enumerate(centers_r):
        _logger.info('For pinhole %i', idx)
        if not cmask[idx]:
            _logger.info('skipping')
            # Fill result with -99
            continue

        # Initial photometric radius
        rad = 3.0
        # Loop to find better photometry radius and background annulus
        irad = rad
        bck = 0.0
        for i in range(phot_niter):
            phot_rad = rad
            # Sky background annulus
            rs1 = rad + back_buff
            rs2 = rs1 + back_width
            _logger.debug('Iter %d, annulus r1=%5.2f r2=%5.2f', i, rs1, rs2)
            bckestim = AnnulusBackgroundEstimator(r1=rs1, r2=rs2)

            # Crop the image to obtain the background
            sl_sky = image_box2d(x0, y0, data.shape, (rs2, rs2))
            raster_sky = data[sl_sky]
            # Logical coordinates
            xx0 = x0 - sl_sky[1].start
            yy0 = y0 - sl_sky[0].start
            # FIXME, perhaps we dont need to crop the image
            try:
                bck = bckestim(raster_sky, xx0, yy0)
                _logger.debug('Iter %d, background %f in '
                              'annulus r1=%5.2f r2=%5.2f',
                              i, bck, rs1, rs2)
            except Exception as error:
                _logger.warning('Error in background estimation %s', error)
                break

            # Radius of the fit
            fit_rad = max(rplot, rad)

            sl = image_box2d(x0, y0, data.shape, (fit_rad, fit_rad))
            part = data[sl]
            # Logical coordinates
            xx0 = x0 - sl[1].start
            yy0 = y0 - sl[0].start

            yy, xx = numpy.mgrid[sl]
            _logger.debug('Iter %d, radial fit', i)

            # Photometry
            dist = numpy.sqrt((xx - x0) ** 2 + (yy - y0) ** 2)
            phot_mask = dist < fit_rad
            r1 = dist[phot_mask]
            part_s = part - bck
            f1 = part_s[phot_mask]
            # Fit radial profile
            model = models.Gaussian1D(amplitude=f1.max(), mean=0, stddev=1.0)
            model.mean.fixed = True  # Mean is always 0.0

            g1d_f = fitter(model, r1, f1, weights=(r1 + 1e-12) ** -1)

            rpeak = g1d_f.amplitude.value
            # sometimes the fit is negative
            rsigma = abs(g1d_f.stddev.value)

            rfwhm = rsigma * FWHM_G

            rad = 2.5 * rfwhm
            _logger.debug('Iter %d, new rad is %f', i, rad)
            if abs(rad - irad) < 1e-3:
                # reached convergence
                _logger.debug('Convergence in iter %d', i)
                break
            else:
                irad = rad
        else:
            _logger.debug('no convergence in photometric radius determination')

        _logger.info('background %6.2f, r1 %7.2f r2 %7.2f', bck, rs1, rs2)
        mm0[idx, 5:5 + 3] = bck, rs1, rs2
        aper_rad = rad
        ca = CircularAperture([(xx0, yy0)], aper_rad)
        m = photutils.aperture_photometry(part_s, ca)
        flux_aper = m['aperture_sum'][0]
        _logger.info('aper rad %f, aper flux %f', aper_rad, flux_aper)
        mm0[idx, 8:8 + 2] = aper_rad, flux_aper

        _logger.info('Radial fit, peak: %f fwhm %f', rpeak, rfwhm)

        try:
            dpeak, dfwhm, smsg = compute_fwhm_enclosed_direct(
                part_s, xx0, yy0, maxrad=fit_rad)
            _logger.info('Enclosed direct, peak: %f fwhm %f', dpeak, dfwhm)
        except Exception as error:
            _logger.warning('Error in compute_fwhm_enclosed_direct %s', error)
            dpeak, dfwhm = -99.0, -99.0

        try:
            eamp, efwhm, epeak, emsg = compute_fwhm_enclosed_grow(
                part_s, xx0, yy0, maxrad=fit_rad)
            _logger.info('Enclosed fit, peak: %f fwhm %f', epeak, efwhm)
        except Exception as error:
            _logger.warning('Error in compute_fwhm_enclosed_grow %s', error)
            eamp, efwhm, epeak, emsg = [-99.0] * 4

        mm0[idx, 10:10 + 6] = epeak, efwhm, dpeak, dfwhm, rpeak, rfwhm

        try:
            res_simple = compute_fwhm_2d_simple(part_s, xx0, yy0)
            _logger.info('Simple, peak: %f fwhm x %f fwhm %f', *res_simple)
            mm0[idx, 16:16 + 3] = res_simple
        except Exception as error:
            _logger.warning('Error in compute_fwhm_2d_simple %s', error)
            mm0[idx, 16:16 + 3] = -99.0

        try:
            res_spline = compute_fwhm_2d_spline(part_s, xx0, yy0)
            _logger.info('Spline, peak: %f fwhm x %f fwhm %f', *res_spline)
            mm0[idx, 19:19 + 3] = res_spline
        except Exception as error:
            _logger.warning('Error in compute_fwhm_2d_spline %s', error)
            mm0[idx, 19:19 + 3] = -99.0

        # Bidimensional fit
        # Fit in a smaller box
        fit2d_rad = int(math.ceil(fit_rad))

        fit2d_half_box = (fit2d_rad, fit2d_rad)
        sl1 = image_box2d(x0, y0, data.shape, fit2d_half_box)

        part1 = data[sl1]
        yy1, xx1 = numpy.mgrid[sl1]

        g2d = models.Gaussian2D(amplitude=rpeak, x_mean=x0, y_mean=y0,
                                x_stddev=1.0, y_stddev=1.0)
        g2d_f = fitter(g2d, xx1, yy1, part1 - bck)

        res_gauss2d = (g2d_f.amplitude.value,
                       g2d_f.x_mean.value + 1,  # FITS coordinates
                       g2d_f.y_mean.value + 1,  # FITS coordinates
                       g2d_f.x_stddev.value * FWHM_G,
                       g2d_f.y_stddev.value * FWHM_G,
                       g2d_f.theta.value
                       )

        _logger.info('Gauss2d, %s', res_gauss2d)
        mm0[idx, 22:22 + 6] = res_gauss2d
        # Moments
        moments_half_box = fit2d_half_box
        res_moments = moments(data, x0, y0, moments_half_box)
        _logger.info('Mxx %f Myy %f Mxy %f e %f pa %f', *res_moments)

        mm0[idx, 28:28 + 5] = res_moments

    # Photometry in coordinates
    # x=centers_r[:,0]
    # y=centers_r[:,1]
    # with radius 2.0 pixels and 4.0 pixels

    apertures = [2.0, 4.0]
    for idx, rad in enumerate(apertures):
        _logger.info('compute photometry with aperture radii %s', rad)
        apertures = CircularAperture(centers_r[:,:2], r=rad)
        faper = photutils.aperture_photometry(data, apertures)
        mm0[:, 33+idx] = faper['aperture_sum']
        _logger.info('done')

    # FITS coordinates
    mm0[:, :4] += 1

    return mm0


def shape_of_slices(tup_of_s):
    return tuple(m.stop - m.start for m in tup_of_s)


def normalize(data):
    b = data.max()
    a = data.min()

    if b != a:
        data_22 = (2 * data - b - a) / (b - a)
    else:
        data_22 = data - b
    return data_22


def normalize_raw(arr):
    """Rescale float image between 0 and 1

    This is an extension of image_as_float of scikit-image
    when the original image was uint16 but later was
    processed and transformed to float32

    Parameters
    ----------
    arr; ndarray

    Returns
    -------
      A ndarray mapped between 0 and 1
    """

    # FIXME: use other limits acording to original arr.dtype
    # This applies only to uint16 images
    # As images were positive, the range is 0,1

    return numpy.clip(arr / 65535.0, 0.0, 1.0)


def char_slit(data, regions, box_increase=3, slit_size_ratio=4.0):

    result = []

    for r in regions:
        _logger.debug('initial region %s', r)
        oshape = shape_of_slices(r)

        ratio = oshape[0] / oshape[1]
        if (slit_size_ratio > 0)  and (ratio < slit_size_ratio):
            _logger.debug("this is not a slit, ratio=%f", ratio)
            continue

        _logger.debug('initial shape %s', oshape)
        _logger.debug('ratio %f', ratio)
        rp = expand_region(r, box_increase, box_increase,
                           start=0, stop=2048)
        _logger.debug('expanded region %r', rp)
        ref = rp[0].start, rp[1].start
        _logger.debug('reference point %r', ref)
    
        datas = data[rp]

        c = ndimage.center_of_mass(datas)
    
        fc = datas.shape[0] // 2
        cc = datas.shape[1] // 2
        _logger.debug("%d %d %d %d", fc, cc, c[0], c[1])

        _peak, fwhm_x, fwhm_y = compute_fwhm_2d_simple(datas, c[1], c[0])

        _logger.debug('x=%f y=%f', c[1] +  ref[1], c[0] +  ref[0])
        _logger.debug('fwhm_x %f fwhm_y %f', fwhm_x, fwhm_y)

        # colrow = ref[1] + cc + 1, ref[0] + fc + 1

        result.append([c[1] +  ref[1] + 1, c[0] +  ref[0] + 1, fwhm_x, fwhm_y])

    return result
