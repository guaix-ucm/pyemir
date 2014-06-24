#
# Copyright 2014 Universidad Complutense de Madrid
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

from __future__ import division

import math

import numpy as np
import scipy.interpolate as itpl
from scipy.interpolate import splrep, splev, sproot
from astropy.modeling import models, fitting
from astropy.modeling import *
from photutils import CircularAnnulus, aperture_circular
from numina.array.mode import mode_half_sample
from numina.array.recenter import (wcs_to_pix_np, img_box, centering_centroid,
        wc_to_pix_1d)


FWHM_G = 2.35482004503

def img_box2d(x, y, shape, box):
    return img_box((y,x), shape, box)

def comp_back_with_annulus(img, xc, yc, rad1, rad2, frac=0.1):
    '''
    center: [x,y], center of first pixel is [0,0]
    '''
    
    ca = CircularAnnulus(rad1, rad2)
    mm = ca.encloses(-0.5 - xc, img.shape[1] - 0.5 - xc, 
                 -0.5 - yc, img.shape[0] - 0.5 - yc, 
                  img.shape[1], img.shape[0])
    
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

# Enclosed gaussian model

class EnclosedGaussian(ParametricModel):
    amplitude = Parameter('amplitude')
    stddev = Parameter('stddev')
    
    def __init__(self, amplitude, stddev, param_dim=1, **constraints):
        super(EnclosedGaussian, self).__init__(
            amplitude=amplitude, stddev=stddev, param_dim=param_dim,
            **constraints)
    
    @staticmethod
    def eval(x, amplitude, stddev):
        return amplitude * (1 - np.exp(-0.5 * (x / stddev)**2))
    
    @staticmethod
    def deriv(x, amplitude, stddev):
        z = (x / stddev)**2
        t = np.exp(-0.5 * z)
        d_amplitude = -t + 1.0
        d_stddev = -amplitude * t * z / stddev
        return [d_amplitude, d_stddev]
    
    @format_input
    def __call__(self, x):
        return self.eval(x, *self.param_sets)

def compute_fwhm_enclosed(imgs, xc, yc, minrad=0.01, maxrad=15.0):

    peak_pix = wcs_to_pix_np((xc, yc))
    peak = imgs[tuple(peak_pix)]
    
    rad = np.logspace(np.log10(minrad), np.log10(maxrad), num=100)
    flux = aperture_circular(imgs, [xc], [yc], rad, method='exact')[:,0]
    
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
    flux = aperture_circular(imgs, [xc], [yc], rad, method='exact')[:,0]
    
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
    flux = aperture_circular(imgs, [xc], [yc], rad, method='exact')[:,0]    
    idx = flux.argmax()
    
    rmodel = rad[:idx+1]
    fmodel = flux[:idx+1]
    fmax = fmodel[-1]
        
    return fit_fwhm_enclosed_grow(fmax, rmodel, fmodel)

def fit_fwhm_enclosed_grow(fmax, rad, flux):
        
    fitter = fitting.NonLinearLSQFitter()
    model1 = EnclosedGaussian(amplitude=fmax, stddev=1.0)
    result = fitter(model1, rad, flux)
        
    amplitude = result.amplitude.value
    sigma = result.stddev.value
    fwhm = sigma * FWHM_G
    peak = amplitude / (2 * math.pi * sigma * sigma)
        
    return amplitude, fwhm, peak, "done"


def compute_fwhm_simple(img, xc, yc):
    
    xpix = wc_to_pix_1d(xc)
    ypix = wc_to_pix_1d(yc)
    
    peak = img[ypix, xpix]
    
    X = range(img.shape[1])
    Y = range(img.shape[0])
    
    res11 = img[ypix, :]
    res22 = img[:, xpix]

    fwhm_x, _codex, _msgx = compute_fwhm_1d(X, res11-0.5 * peak, xc, xpix)
    fwhm_y, _codey, _msgy = compute_fwhm_1d(Y, res22-0.5 * peak, yc, ypix)
    
    return peak, fwhm_x, fwhm_y


def compute_fwhm_spline2d_fit(img, xc, yc):

    xpix = wc_to_pix_1d(xc)
    ypix = wc_to_pix_1d(yc)
    # The image is already cropped
    Y = np.arange(img.shape[1])
    X = np.arange(img.shape[0])
    bb = itpl.RectBivariateSpline(X, Y, img)
    # We assume that the peak is in the center...
    peak = bb(xc, yc)[0,0]
    
    U = X
    V = bb.ev(U, [yc for _ in U]) - 0.5 * peak
    fwhm_x, codex, msgx = compute_fwhm_1d(U, V, yc, ypix)

    U = Y
    V = bb.ev([xc for _ in U], U) - 0.5 * peak
    fwhm_y, codey, msgy = compute_fwhm_1d(U, V, xc, xpix)

    return peak, fwhm_x, fwhm_y



def _fwhm_side_lineal(uu, vv):
    '''Compute r12 using linear interpolation.'''
    res1, = np.nonzero(vv < 0)
    if len(res1) == 0:
        return 0, 1 # error, no negative value
    else:
        # first value
        i2 = res1[0]
        i1 = i2 -1
        dx = uu[i2] - uu[i1]
        dy = vv[i2] - vv[i1]
        r12 = uu[i1] - vv[i1] * dx / dy
        return r12, 0

def compute_fwhm_1d(uu, vv, uc, upix):

    _fwhm_side = _fwhm_side_lineal
    
    # Find half peak radius on the rigth
    r12p, errorp = _fwhm_side(uu[upix:], vv[upix:])
    
    # Find half peak radius on the left
    r12m, errorm = _fwhm_side(uu[upix::-1], vv[upix::-1])
    
    if errorm == 1:
        if errorp == 1:
            fwhm = -99 # No way
            msg = 'Failed to compute FWHM'
            code = 2
        else:
            fwhm = 2 * (r12p - uc)
            code = 1
            msg = 'FWHM computed from rigth zero'
    else:
        if errorp == 1:
            fwhm = 2 * (uc - r12m)
            msg = 'FWHM computed from left zero'
            code = 1
        else:
            msg = 'FWHM computed from left and rigth zero'
            code = 0
            fwhm = r12p - r12m

    
    return fwhm, code, msg

def extent(sl):
    result = [sl[1].start-0.5, sl[1].stop-0.5, sl[0].start-0.5, sl[0].stop-0.5]
    return result

def moments(data, x0, y0, half_box):
    
    sl1 = img_box2d(x0, y0, data.shape, half_box)
    
    part1 = data[sl1]
    norm = part1.sum()
    
    Y1, X1 = np.mgrid[sl1]
    
    xo = X1 - x0
    yo = Y1 - y0
    
    Mxx = (xo**2 * part1).sum() / norm
    Myy = (yo**2 * part1).sum() / norm
    Mxy = (xo * yo * part1).sum() / norm
    e = math.sqrt((Mxx - Myy)**2 + (2 * Mxy)** 2) / (Mxx + Myy)
    pa = 0.5 * math.atan(2 * Mxy / (Mxx - Myy))
    
    return Mxx, Myy, Mxy, e, pa


def rim(data, xinit, yinit, 
        recenter_half_box=5, 
        recenter_nloop=10, 
        recenter_maxdist=10.0, buff=3, width=5, niter=10, rplot=8):

    fine_recentering = False
    sigma0 = 1.0
    rad =  3 * sigma0 * FWHM_G
    box = recenter_half_box
    recenter_half_box = (box, box)
    plot_half_box = (3*box, 3*box)
    
    print 'C initial', xinit, yinit
    x0, y0, back_, msg = centering_centroid(data, xinit, yinit, 
     box=recenter_half_box,
     maxdist=recenter_maxdist, nloop=recenter_nloop)
    print 'C final', x0, y0
    
    if fine_recentering:
        print 'Fine recentering'
        print 'C initial', x0, y0
        x1, y1, back_, msg = centering_centroid(data, x0, y0, box=(1,1), 
                maxdist=2*math.sqrt(2), nloop=1)
        print 'C final', x1, y1
    
    sl = img_box2d(x0, y0, data.shape, plot_half_box)
    part = data[sl]
    
    
    xx0 = x0 - sl[1].start
    yy0 = y0 - sl[0].start
    
    Y, X = np.mgrid[sl]
    
    # Photometry
    D = np.sqrt((X-x0)**2 + (Y-y0)**2)
    m = D < rplot
    r1 = D[m]
    
    fitter = fitting.NonLinearLSQFitter()
    model = models.Gaussian1D(amplitude=1.0, mean=0, stddev=1.0)
    model.mean.fixed = True # Mean is always 0.0
        
    irad = rad
    for i in range(niter):
        rs1 = rad + buff
        rs2 = rs1 + width
        
        bckestim = AnnulusBackgroundEstimator(r1=rs1, r2=rs2)
        bck = bckestim(part, xx0, yy0)
        part_s = part - bck

        flux_aper = aperture_circular(part_s, [xx0], [yy0], rad, method='exact')
    
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
            print 'convergence in iter %d' % (i+1)
            break
        else:
            irad = rad
    else:
        print 'no convergence in photometric radius determination'


    print 'P, aper rad', rad, 'flux_aper', flux_aper[0]
    print 'P, annulus background:', bck, 'radii', rs1, rs2
    eamp, efwhm, epeak, emsg = compute_fwhm_enclosed_grow(part_s, xx0, yy0, maxrad=rs1)
    print 'Enclosed fit, peak:', epeak, 'fwhm', efwhm
    print 'Radial fit, peak:', rpeak, 'fwhm', rfwhm
    print 'Direct enclosed, peak:', dpeak, 'dfwhm', dfwhm
    
    lpeak, fwhm_x, fwhm_y = compute_fwhm_simple(part_s, xx0, yy0)
    print 'Simple, peak:', lpeak, 'fwhm x', fwhm_x, 'fwhm y', fwhm_y
   
    # Fit in a smaller box
    fit2d_rad = int(math.ceil(0.5 * rad))
    
    fit2d_half_box = (fit2d_rad, fit2d_rad)
    sl1 = img_box2d(x0, y0, data.shape, fit2d_half_box)
    part1 = data[sl1]
    Y1, X1 = np.mgrid[sl1]
    g2d = models.Gaussian2D(amplitude=rpeak, x_mean=x0, y_mean=y0, x_stddev=1.0, y_stddev=1.0)
    g2d_f = fitter(g2d, X1, Y1, part1 - bck)
    print 'Gauss2D fit'
    print g2d_f
    
    moments_half_box = fit2d_half_box
    Mxx, Myy, Mxy, e, pa = moments(data, x0, y0, moments_half_box)

    print Mxx, Myy, Mxy, e, pa
