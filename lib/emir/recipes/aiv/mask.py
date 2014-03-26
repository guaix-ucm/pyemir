#
# Copyright 2013-2014 Universidad Complutense de Madrid
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

import logging
import math

import numpy

from astropy.io import fits
import photutils
from astropy.modeling import models, fitting
from scipy.spatial import distance
import scipy.interpolate as itpl
import scipy.optimize as opz

import scipy.ndimage.filters as filters
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import scipy.stats
from scipy.spatial import cKDTree as KDTree
import scipy.spatial.distance
import scipy.ndimage.measurements
from scipy.ndimage.interpolation import map_coordinates
from scipy.interpolate import UnivariateSpline
import scipy.interpolate as itpl
import scipy.optimize as opz


from numina.core import RecipeError
from numina.core import BaseRecipe, RecipeRequirements, DataFrame
from numina.core import Requirement, Product, DataProductRequirement, Parameter
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement

from numina.array.combine import median, mean
from numina import __version__
from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import FlatFieldCorrector, SkyCorrector
from numina.flow import SerialFlow
from numina.flow.processing import DivideByExposure

from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import FrameDataProduct, MasterIntensityFlat
from emir.dataproducts import DarkCurrentValue, CoordinateList2DType
from emir.dataproducts import CentroidsTableType, ArrayType

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"
            
GAUSS_FWHM_FACTOR = 2.354800

def img_box(center, shape, box):

    def slice_create(c, s, b):
        cc = int(math.floor(c + 0.5))
        l = max(0, cc - b)
        h = min(s, cc + b +1)
        return slice(l, h, None)

    return tuple(slice_create(*args) for args in zip(center, shape, box))

# returns y,x
def _centering_centroid_loop(data, center, box):
    sl = img_box(center, data.shape, box)
    raster = data[sl]
    mm = raster.mean()
    
    threshold = mm
    _logger.debug('Threshold is %s', threshold)
    
    mask = raster >= threshold
    if not numpy.any(mask):
        _logger.warning('No points to compute centroid, threshold too high')
        return center
        
    rr = numpy.where(mask, raster, 0)

    r_std = rr.std()
    r_mean = rr.mean
    if r_st > 0:
        snr = r_mean / r_std
        _logger.debug('SNR is %f', snr)
    
    fi, ci = numpy.indices(raster.shape)
    
    norm = rr.sum()
    print norm
    fm = (rr * fi).sum() / norm
    cm = (rr * ci).sum() / norm
    
    return fm + sl[0].start, cm + sl[1].start
    
# returns y,x
def centering_centroid(data, center, box, nloop=10):
    toldist = 1e-3
    ncenter = center
    for i in range(nloop):
        center = _centering_centroid_loop(data, ncenter, box)
        # check convergence
        dst = distance.euclidean(ncenter, center)
        if dst < toldist:
            return ncenter, 'converged in iteration %i' % i
        else:
            ncenter = center
    
    return ncenter, 'not converged in %i iterations' % nloop

# returns y,x
def compute_fwhm(img, center):
    X = numpy.arange(img.shape[0])
    Y = numpy.arange(img.shape[1])

    bb = itpl.RectBivariateSpline(X, Y, img)
    # We assume that the peak is in the center...
    peak = bb(*center)[0,0]

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
        rV = V[cp-1::-1]
        rU = U[cp-1::-1]
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

    return peak, center[0], center[1], peak, fwhm_x, fwhm_y

# returns x,y
def compute_fwhm_global(data, center, box):
    sl = img_box(center, data.shape, box)
    raster = data[sl]
    newc = center[0] - sl[0].start, center[1] - sl[1].start
    res = compute_fwhm(raster, newc)
    return res[1]+sl[1].start, res[0]+sl[0].start, res[2], res[3], res[4]
    #return res[0]+sl[0].start, res[1]+sl[1].start, res[2], res[3], res[4]

# returns x,y
def gauss_model(data, center_r):
    sl = img_box(center_r, data.shape, box=(4,4))
    raster = data[sl]
    
    new_c = center_r[0] - sl[0].start, center_r[1] - sl[1].start

    # background
    mm = raster - raster.min()

    yi, xi = numpy.indices(raster.shape)

    g = models.Gaussian2D(amplitude=1.2, x_mean=new_c[1], y_mean=new_c[0], x_stddev=1.0, y_stddev=1.0)
    f1 = fitting.NonLinearLSQFitter()
    t = f1(g, xi, yi, mm)
    mm = t.x_mean.value + sl[1].start, t.y_mean.value+ sl[0].start, t.amplitude.value, t.x_stddev.value, t.y_stddev.value
    return mm

def pinhole_char(data, ncenters, box=4):

    ibox = (box, box)
    
    # convert FITS x, y coordinates (pixel center 1)
    # to python/scipy/astropy (pixel center 0 and x,y -> y,x)
    
    centers_py = numpy.fliplr(ncenters[:,0:2]) - 1
    
    # recentered values
    centers_r = numpy.empty_like(centers_py)
    
    for idx, c in enumerate(centers_py):
        center, _msg = centering_centroid(data, c, box=ibox)
        centers_r[idx] = center
        
    mm1 = numpy.empty((centers_r.shape[0], 5))
    
    # compute the FWHM without fitting
    for idx, center in enumerate(centers_r):
        res = compute_fwhm_global(data, center, box=ibox)
        mm1[idx] = res
    
    mm2 = numpy.empty((centers_r.shape[0], 5))
    # compute the FWHM fitting a Gaussian
    for idx, center in enumerate(centers_r):
        mm2[idx] = gauss_model(data, center)
    
    # Photometry in coordinates
    # x=centers_r[:,1]
    # y=centers_r[:,0]
    # with radius 2.0 pixels and 4.0 pixels
    mm3 = photutils.aperture_circular(data, centers_r[:,1], centers_r[:,0], [2., 4])
    
    # Convert coordinates to FITS
    mm1[:,0:2] += 1
    mm2[:,0:2] += 1
    
    return mm1, mm2, mm3

class TestPinholeRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')
    master_sky = DataProductRequirement(MasterIntensityFlat, 'Master Sky calibration')
    pinhole_nominal_positions = Requirement(CoordinateList2DType, 'Nominal positions of the pinholes')
    shift_coordinates = Parameter(True, 'Use header information to shift the pinhole positions from (0,0) to X_DTU, Y_DTU')
    box_half_size = Parameter(4, 'Half of the search box size in pixels')

class TestPinholeRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)
    fwhm = Product(ArrayType)
    gauss = Product(ArrayType)
    phot = Product(ArrayType)
    

@define_requirements(TestPinholeRecipeRequirements)
@define_result(TestPinholeRecipeResult)
class TestPinholeRecipe(BaseRecipe):

    def __init__(self):
        super(TestPinholeRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting pinhole processing')

        # Loading calibrations
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            mbias = hdul[0].data
            bias_corrector = BiasCorrector(mbias)
            #bias_corrector = IdNode()
            
        exposure_corrector = DivideByExposure(factor=1.0)

        with rinput.master_dark.open() as mdark_hdul:
            _logger.info('loading dark')
            mdark = mdark_hdul[0].data
            dark_corrector = DarkCorrector(mdark)

        with rinput.master_flat.open() as mflat_hdul:
            _logger.info('loading intensity flat')
            mflat = mflat_hdul[0].data
            flat_corrector = FlatFieldCorrector(mflat)

        with rinput.master_sky.open() as msky_hdul:
            _logger.info('loading sky')
            msky = msky_hdul[0].data
            sky_corrector = SkyCorrector(msky)

        flow = SerialFlow([bias_corrector, exposure_corrector, 
                dark_corrector, flat_corrector, sky_corrector])

        odata = []
        cdata = []        
        try:
            _logger.info('processing input frames')
            for frame in rinput.obresult.frames:                
                hdulist = frame.open()
                final = flow(hdulist)
                cdata.append(final)
                odata.append(hdulist)

            if len(cdata) > 1:
                _logger.info('stacking %d images using median', len(cdata))
                data = median([d['primary'].data for d in cdata], dtype='float32')
                hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header.copy())
            else:
                _logger.info('skip stacking 1 image')
                img = cdata[0]
                hdu = img[0].copy()

        finally:
            _logger.debug('closing images')
            for hdulist in odata:
                hdulist.close()
                
            del cdata
        
            
        _logger.debug('update result header')
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdulist = fits.HDUList([hdu])

        _logger.debug('finding pinholes')
        
        if rinput.shift_coordinates:
            #get things from header
            _logger.info('getting DTU position from header')
            try:
                xdtu = hdr.get('X_DTU')
                ydtu = hdr.get('Y_DTU')
                zdtu = hdr.get('Z_DTU')
            except KeyError as error:
                _logger.error(error)
                raise RecipeError
            _logger.info('X_DTU=%6.2f Y_DTU=%6.2f Z_DTU=%6.2f', xdtu, ydtu, zdtu)
            # transform coordinates
            _logger.info('transform pinhole coordinates from reference (0,0)')
            xfac = xdtu * 0.055
            yfac = -ydtu * 0.055
        
            vec = numpy.array([yfac, xfac])
            _logger.info('shift is %s', vec)
            ncenters = rinput.pinhole_nominal_positions + vec
        else:
            _logger.info('using pinhole coordinates as they are')
            ncenters = rinput.pinhole_nominal_positions        
        
        
        fwhm, gauss, phot = pinhole_char(hdul[0].data, ncenters, box=rinput.box_half_size)
        
        result = TestPinholeRecipeResult(frame=hdul, fwhm=fwhm, gauss=gauss, phot=phot)
        return result
        
        
