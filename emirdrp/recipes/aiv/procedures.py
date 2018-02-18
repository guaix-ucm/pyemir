#
# Copyright 2014-2015 Universidad Complutense de Madrid
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

'''AIV Recipes for EMIR'''

from __future__ import division, print_function

import math

import numpy as np
from scipy.interpolate import splrep, splev, sproot
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils.geometry.circular_overlap import circular_overlap_grid

from astropy.modeling import (fitting, models)

from numina.array.mode import mode_half_sample
from numina.array.recenter import centering_centroid
from numina.array.utils import wcs_to_pix_np
from numina.array.fwhm import compute_fwhm_1d_simple
from numina.array.utils import image_box2d
from numina.modeling import EnclosedGaussian
from numina.constants import FWHM_G


def encloses_annulus(x_min, x_max, y_min, y_max, nx, ny, r_in, r_out):
    '''Encloses function backported from old photutils'''

    gout = circular_overlap_grid(x_min, x_max, y_min, y_max, nx, ny, r_out, 1, 1)
    gin = circular_overlap_grid(x_min, x_max, y_min, y_max, nx, ny, r_in, 1, 1)
    return gout - gin


def comp_back_with_annulus(img, xc, yc, r_in, r_out, frac=0.1):
    '''
    center: [x,y], center of first pixel is [0,0]
    '''

    x_min = -0.5 - xc
    x_max = img.shape[1] - 0.5 - xc
    y_min = -0.5 - yc
    y_max = img.shape[1] - 0.5 - yc
    mm = encloses_annulus(x_min, x_max, y_min, y_max,
                          img.shape[1], img.shape[0],
                          r_in, r_out)

    valid = mm > frac
    rr = img[valid]

    if rr.size == 0:
        raise ValueError("Not enough points to compute background")

    # mode?
    bck = mode_half_sample(rr)

    return bck, mm


class ConstantBackgroundEstimator(object):
    def __init__(self, const):
        self.const = const

    def __call__(self, a, x, y):
        return self.const


class AnnulusBackgroundEstimator(object):
    def __init__(self, r1, r2):
        self.r1 = r1
        self.r2 = r2

    def __call__(self, a, x, y):
        bck, _ = comp_back_with_annulus(a, x, y, self.r1, self.r2)
        return bck


def compute_fwhm_enclosed(imgs, xc, yc, minrad=0.01, maxrad=15.0):

    peak_pix = wcs_to_pix_np((xc, yc))
    peak = imgs[tuple(peak_pix)]

    rad = np.logspace(np.log10(minrad), np.log10(maxrad), num=100)
    flux = np.zeros_like(rad)
    positions = [(xc, yc)]
    for idr, r in enumerate(rad):
        ca = CircularAperture(positions, r)
        m = aperture_photometry(imgs, ca)
        flux[idr] = m['aperture_sum'][0]

    idx = flux.argmax()

    rmodel = rad[:idx+1]
    fmodel = flux[:idx+1]
    fmax = fmodel[-1]

    res_d = fit_fwhm_enclosed_direct(peak, rad, flux)

    res_g = fit_fwhm_enclosed_grow(fmax, rmodel, fmodel)

    return res_d, res_g


def compute_fwhm_enclosed_direct(imgs, xc, yc, minrad=0.01, maxrad=15.0):

    peak_pix = wcs_to_pix_np((xc, yc))
    peak = imgs[tuple(peak_pix)]

    rad = np.logspace(np.log10(minrad), np.log10(maxrad), num=100)
    flux = np.zeros_like(rad)
    positions = [(xc, yc)]
    for idr, r in enumerate(rad):
        ca = CircularAperture(positions, r)
        m = aperture_photometry(imgs, ca)
        flux[idr] = m['aperture_sum'][0]

    return fit_fwhm_enclosed_direct(peak, rad, flux)


def fit_fwhm_enclosed_direct(peak, rad, flux):

    # We use splines to interpolate and derivate
    spl = splrep(rad, flux)
    # First derivative
    vald1 = splev(rad, spl, der=1)

    splinter = splrep(rad, vald1 - math.pi * peak * rad)

    roots = sproot(splinter)
    nroots = len(roots)
    if peak < 0:
        msg = "The method doesn't converge, peak is negative"
        fwhm = -99
    else:
        if nroots == 0:
            msg = "The method doesn't converge, no roots"
            fwhm = -99
        elif nroots == 1:
            r12 = roots[0]
            fwhm = 2 * r12
            msg = "The method converges, one root"
        else:
            msg = "The method doesn't converge, multiple roots"
            r12 = roots[0]
            fwhm = 2 * r12

    return peak, fwhm, msg


def compute_fwhm_enclosed_grow(imgs, xc, yc, minrad=0.01, maxrad=15.0):

    rad = np.logspace(np.log10(minrad), np.log10(maxrad), num=100)
    flux = np.zeros_like(rad)
    positions = [(xc, yc)]
    for idr, r in enumerate(rad):
        ca = CircularAperture(positions, r)
        m = aperture_photometry(imgs, ca)
        flux[idr] = m['aperture_sum'][0]
    idx = flux.argmax()
    rmodel = rad[:idx+1]
    fmodel = flux[:idx+1]
    fmax = fmodel[-1]

    return fit_fwhm_enclosed_grow(fmax, rmodel, fmodel)


def fit_fwhm_enclosed_grow(fmax, rad, flux):

    fitter = fitting.LevMarLSQFitter()
    model1 = EnclosedGaussian(amplitude=fmax, stddev=1.0)
    result = fitter(model1, rad, flux)

    amplitude = result.amplitude.value
    sigma = result.stddev.value
    fwhm = sigma * FWHM_G
    peak = amplitude / (2 * math.pi * sigma * sigma)

    return amplitude, fwhm, peak, "done"


def moments(data, x0, y0, half_box):

    sl1 = image_box2d(x0, y0, data.shape, half_box)

    part1 = data[sl1]
    norm = part1.sum()

    Y1, X1 = np.mgrid[sl1]

    xo = X1 - x0
    yo = Y1 - y0

    Mxx = (xo**2 * part1).sum() / norm
    Myy = (yo**2 * part1).sum() / norm
    Mxy = (xo * yo * part1).sum() / norm
    e = math.sqrt((Mxx - Myy)**2 + (2 * Mxy)**2) / (Mxx + Myy)
    pa = 0.5 * math.atan(2 * Mxy / (Mxx - Myy))

    return Mxx, Myy, Mxy, e, pa


def rim(data, xinit, yinit,
        recenter_half_box=5,
        recenter_nloop=10,
        recenter_maxdist=10.0, buff=3, width=5, niter=10, rplot=8):

    fine_recentering = False
    sigma0 = 1.0
    rad = 3 * sigma0 * FWHM_G
    box = recenter_half_box
    recenter_half_box = (box, box)
    plot_half_box = (3*box, 3*box)

    print('C initial', xinit, yinit)
    x0, y0, _1, _2, _3 = centering_centroid(data, xinit, yinit,
                                            box=recenter_half_box,
                                            maxdist=recenter_maxdist,
                                            nloop=recenter_nloop
                                            )
    print('C final', x0, y0)

    if fine_recentering:
        print('Fine recentering')
        print('C initial', x0, y0)
        x1, y1, _back, _status, _msg = centering_centroid(
            data,
            x0,
            y0,
            box=(1, 1),
            maxdist=2*math.sqrt(2),
            nloop=1
            )
        print('C final', x1, y1)

    sl = image_box2d(x0, y0, data.shape, plot_half_box)
    part = data[sl]

    xx0 = x0 - sl[1].start
    yy0 = y0 - sl[0].start

    Y, X = np.mgrid[sl]

    # Photometry
    D = np.sqrt((X-x0)**2 + (Y-y0)**2)
    m = D < rplot
    r1 = D[m]

    fitter = fitting.LevMarLSQFitter()
    model = models.Gaussian1D(amplitude=1.0, mean=0, stddev=1.0)
    model.mean.fixed = True  # Mean is always 0.0

    irad = rad
    for i in range(niter):
        rs1 = rad + buff
        rs2 = rs1 + width

        bckestim = AnnulusBackgroundEstimator(r1=rs1, r2=rs2)
        bck = bckestim(part, xx0, yy0)
        part_s = part - bck

        ca = CircularAperture([(xx0, yy0)], rad)
        m = aperture_photometry(part_s, ca)
        flux_aper = m['aperture_sum'][0]

        f1 = part_s[m]
        g1d_f = fitter(model, r1, f1, weights=(r1+1e-12)**-1)

        rpeak = g1d_f.amplitude.value
        # sometimes the fit is negative
        rsigma = abs(g1d_f.stddev.value)

        rfwhm = rsigma * FWHM_G

        dpeak, dfwhm, smsg = compute_fwhm_enclosed_direct(part_s, xx0, yy0)

        rad = 3 * dfwhm
        if abs(rad-irad) < 1e-3:
            # reached convergence
            print('convergence in iter %d' % (i+1))
            break
        else:
            irad = rad
    else:
        print('no convergence in photometric radius determination')

    print('P, aper rad', rad, 'flux_aper', flux_aper[0])
    print('P, annulus background:', bck, 'radii', rs1, rs2)
    eamp, efwhm, epeak, emsg = compute_fwhm_enclosed_grow(
        part_s, xx0, yy0, maxrad=rs1
        )
    print('Enclosed fit, peak:', epeak, 'fwhm', efwhm)
    print('Radial fit, peak:', rpeak, 'fwhm', rfwhm)
    print('Direct enclosed, peak:', dpeak, 'dfwhm', dfwhm)

    lpeak, fwhm_x, fwhm_y = compute_fwhm_1d_simple(part_s, xx0, yy0)
    print('Simple, peak:', lpeak, 'fwhm x', fwhm_x, 'fwhm y', fwhm_y)

    # Fit in a smaller box
    fit2d_rad = int(math.ceil(0.5 * rad))

    fit2d_half_box = (fit2d_rad, fit2d_rad)
    sl1 = image_box2d(x0, y0, data.shape, fit2d_half_box)
    part1 = data[sl1]
    Y1, X1 = np.mgrid[sl1]
    g2d = models.Gaussian2D(amplitude=rpeak, x_mean=x0, y_mean=y0,
                            x_stddev=1.0, y_stddev=1.0)
    g2d_f = fitter(g2d, X1, Y1, part1 - bck)
    print('Gauss2D fit')
    print(g2d_f)

    moments_half_box = fit2d_half_box
    Mxx, Myy, Mxy, e, pa = moments(data, x0, y0, moments_half_box)

    print(Mxx, Myy, Mxy, e, pa)
