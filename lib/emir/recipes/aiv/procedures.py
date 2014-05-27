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
#import scipy.interpolate as itpl
#import scipy.optimize as opz
from scipy.interpolate import splrep, splev, sproot

#from astropy.modeling import models, fitting
#from astropy.io import fits
#import photutils

from photutils import CircularAnnulus, aperture_circular

#from scipy.interpolate import splrep, splev, sproot
from numina.array.mode import mode_half_sample
from numina.array.recenter import wcs_to_pix_np, img_box, centering_centroid
            
#GAUSS_FWHM_FACTOR = 2.354800

FWHM_G = 2.35482004503

def comp_back_with_annulus(img, xc, yc, rad1, rad2, frac=0.1):
    '''
    center: [x,y], center of first pixel is [0,0]
    '''
    
    ca = CircularAnnulus(rad1, rad2)
    mm = ca.encloses(-0.5 - xc, img.shape[1] - 0.5 - xc, 
                 -0.5 - yc, img.shape[0] - 0.5 - yc, 
                  img.shape[0], img.shape[1])
    
    valid = mm > frac
    rr = img[valid]
    
    if rr.size == 0:
        raise ValueError # Bad
    
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

#def enclosed_fwhm(data, x0, y0, bckestim):
#    
#    bck = bckestim(data, x0, y0)
#    #print 'background=', bck
#    
#    peak, fwhm, msg = compute_fwhm_enclosed(data - bck, x0, y0)
#    
#    return x0, y0, peak, bck, fwhm, msg

def compute_fwhm_enclosed(data, x0, y0, maxrad=15.0):
    
    res = data
    peak_pix = wcs_to_pix_np((x0, y0))
    peak = res[tuple(peak_pix)]
    
    min_log_rad = -2
    
    rad = np.logspace(min_log_rad, np.log10(maxrad), num=100)
    flux = aperture_circular(res, [x0], [y0], rad, method='exact')[:,0]

    # We use splines to interpolate and derivate
    spl = splrep(rad, flux)
    vald1 = splev(rad, spl, der=1)
    splinter = splrep(rad, vald1 - math.pi * peak * rad)
    #valdiff = splev(rad, splinter)
    # Evaluate the spline representation
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
            msg = "The method converges"
        else:
            msg = "The method doesn't converge, multiple roots"
            r12 = roots[0]
            fwhm = 2 * r12
    
    return peak, fwhm, msg

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

def compute_fwhm_no(img, xc, yc, bckestim):
    # background estimation
    bck = bckestim(img, xc, yc)
    
    img_b = img - bck
    
    ypix, xpix  = wcs_to_pix_np([xc, yc])
    
    peak = img_b[ypix, xpix]
    
    xitp1 = range(img.shape[1])
    res11 = img_b[ypix, :]
    
    yitp1 = range(img.shape[0])
    res22 = img_b[:, xpix]

    fwhm_x, codex, msgx = compute_fwhm_1d(xitp1, res11-0.5 * peak, xc, xpix)
    fwhm_y, codey, msgy = compute_fwhm_1d(yitp1, res22-0.5 * peak, yc, ypix)
    
    return xc, yc, peak, bck, fwhm_x, fwhm_y
    
from astropy.modeling import *

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

from astropy.modeling import models, fitting

def compute_fwhm_enclosed2(res, x0, y0, maxrad=15.0, extra=False):
    
    log_min_rad = -2
    
    rad = np.logspace(log_min_rad, np.log10(maxrad), num=100)
    flux = aperture_circular(res, [x0], [y0], rad)[:,0]
    
    idx = flux.argmax()
    #print rad[idx], flux[idx]
    #print 'maxflux is', flux[idx], 'at radius', rad[idx], 'index', idx

    #nflux = np.where(rad > rad[idx], flux[idx], flux)
    
    rmodel = rad[:idx+1]
    fmodel = flux[:idx+1]
    
    fitter = fitting.NonLinearLSQFitter()
    model1 = EnclosedGaussian(amplitude=flux[-1], stddev=1.0)
    result1 = fitter(model1, rmodel, fmodel)
    
    #model2 = EnclosedGaussian(amplitude=flux[idx], stddev=1.0)
    #model2.amplitude.fixed = True
    #result2 = fitter(model2, rmodel, fmodel)
    
    amplitude = result1.amplitude.value
    sigma = result1.stddev.value
    fwhm = sigma * FWHM_G
    peak  = amplitude / (2 * math.pi * sigma * sigma)
    
    #plt.plot(rmodel, fmodel, 'ko-')
    #plt.plot(rmodel, result1(rmodel), 'r-')
    ##plt.plot(rmodel, result2(rmodel), 'g-')
    #plt.show()
    
    
    return amplitude, fwhm, peak, "done"

def raster(data, x, y, half_box):
    sl = img_box((y,x), data.shape, half_box)
    extent = [sl[1].start-0.5, sl[1].stop-0.5, sl[0].stop-0.5, sl[0].start-0.5]
    return data[sl], extent

def rim(data, xinit, yinit, 
        recenter_half_box=5, 
        recenter_nloop=10, 
        recenter_maxdist=10.0, buff=3, width=5, niter=3, rplot=8):

    fine_recentering = False
    sigma0 = 2
    rad =  sigma0 * FWHM_G
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
    
    
    part, extent = raster(data, x0, y0, plot_half_box)

    
    # logical coordinates
    xx0 = x0 - (extent[0] + 0.5)
    yy0 = y0 - (extent[3] + 0.5)
    if fine_recentering:
        xx1 = x1 - (extent[0] + 0.5)
        yy1 = y1 - (extent[3] + 0.5)

    #plt.imshow(part, interpolation='nearest')
    #circ = plt.Circle((xx0, yy0), radius=1.0, color='r', alpha=0.5)
    #plt.gca().add_patch(circ)
    #if fine_recentering:
    #    circ = plt.Circle((xx1, yy1), radius=1.0, color='g', alpha=0.5)
    #    plt.gca().add_patch(circ)
    #plt.show()
    
    shape = part.shape
    Y, X = np.mgrid[0:shape[0], 0:shape[1]]
    D = np.sqrt((X-xx0)**2 + (Y-yy0)**2)
    m = D < rplot
    r1 = D[m]
    
    for _ in range(niter):
        rs1 = rad + buff
        rs2 = rs1 + width
        
        bckestim = AnnulusBackgroundEstimator(r1=rs1, r2=rs2)
        bck = bckestim(part, xx0, yy0)
        part_s = part - bck
    
        flux_aper = aperture_circular(part_s, [xx0], [yy0], rad, method='exact')
    
        f1 = part_s[m]
        model = models.Gaussian1D(amplitude=1.0, mean=0, stddev=1.0)
        model.mean.fixed = True # Mean is always 0.0
    
        fitter = fitting.NonLinearLSQFitter()
    
        g = fitter(model, r1, f1, weights=(r1+1e-12)**-1)

        peak = g.amplitude.value
        sigma = g.stddev.value
        fwhm = sigma * FWHM_G

        dpeak, dfwhm, smsg = compute_fwhm_enclosed(part_s, xx0, yy0)
        eamp, efwhm, epeak, emsg = compute_fwhm_enclosed2(part_s, xx0, yy0, maxrad=rs1)
        rad = 3 * dfwhm

#    plt.plot(r1, f1, 'ko');
#    r = np.linspace(0, rplot)
#    plt.plot(r, g(r), 'r-', lw=2, label='Gaussian')
#    plt.show()

    #print 'L', compute_fwhm(part, (xx0, yy0), bckestim)

    print 'P, aper rad', rad, 'flux_aper', flux_aper[0]
    print 'P, annulus background:', bck, 'radii', rs1, rs2
    print 'Enclosed fit, peak:', epeak, 'fwhm', efwhm
    print 'Radial fit, peak:', peak, 'fwhm', fwhm
    print 'Direct enclosed, peak:', dpeak, 'dfwhm', dfwhm
    
    xc, yc, peak, bck, fwhm_x, fwhm_y = compute_fwhm_no(part, xx0, yy0, bckestim)
    print 'L, peak:', peak, 'fwhm x', fwhm_x, 'fwhm y', fwhm_y

