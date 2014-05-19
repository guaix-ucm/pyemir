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
import scipy.interpolate as itpl
import scipy.optimize as opz
from astropy.modeling import models, fitting
from astropy.io import fits
import photutils

from scipy.interpolate import splrep, splev, sproot
from numina.array.mode import mode_half_sample
from numina.array.recenter import wcs_to_pix_np, img_box, centering_centroid_xy

from numina import __version__
from numina.core import BaseRecipe, RecipeRequirements, RecipeError
from numina.core import Requirement, Product, DataProductRequirement, Parameter
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement
from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import FlatFieldCorrector, SkyCorrector
from numina.flow import SerialFlow
from numina.flow.processing import DivideByExposure
from numina.flow.node import IdNode
from numina.array import combine

from emir.core import RecipeResult
from emir.core import EMIR_BIAS_MODES
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import FrameDataProduct, MasterIntensityFlat
from emir.dataproducts import CoordinateList2DType
from emir.dataproducts import ArrayType
from emir.core import gather_info_frames, gather_info_dframe

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"
            
#GAUSS_FWHM_FACTOR = 2.354800

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

    return center[0], center[1], peak, fwhm_x, fwhm_y


# Background in an annulus, mode is HSM

from photutils import CircularAnnulus

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

def compute_fwhm_enclosed(data, x0, y0, back, extra=False):
    
    res = data - back
    #print 'center is', x0, y0
    peak_pix = wcs_to_pix_np((x0, y0))
    #print 'peak pixel is',peak_pix

    peak = res[tuple(peak_pix)]
    
    
    minrad = 0.01
    maxrad = 15.0
    rad = numpy.logspace(-2, numpy.log10(maxrad), num=100)
    flux = photutils.aperture_circular(res, [x0], [y0], rad, method='exact')[:,0]

    # We use splines to interpolate and derivate
    spl = splrep(rad, flux)
    # Evaluate the spline representation
    val = splev(rad, spl)
    # Evaluate the derivative
    vald1 = splev(rad, spl, der=1)

    spl = splrep(rad, flux)
    val = splev(rad, spl)
    vald1 = splev(rad, spl, der=1)
    splinter = splrep(rad, vald1 - math.pi * peak * rad)
    valdiff = splev(rad, splinter)
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

    if extra:
        return peak, fwhm, msg, {'rad': rad, 'flux': flux, 'deriv': vald1}
    else:
        return peak, fwhm, msg

def compute_fwhm_global_no(data, center, box):
    sl = img_box(center, data.shape, box)
    
    raster = data[sl]
    
    bck_estim = AnnulusBackgroundEstimator(r1=5.0, r2=13.0)
    
    y0 = center[0] - sl[0].start 
    x0 = center[1] - sl[1].start
    
    x, y, peak, bck, fwhm_x, fwhm_y = compute_fwhm_no(raster, x0, y0, bck_estim)
    
    res0 = x+sl[1].start, y+sl[0].start, peak, bck, fwhm_x, fwhm_y
    
    x, y, peak, bck, fwhm, _msg  = enclosed_fwhm(data, x0, y0, bck_estim)
    
    res1 = x+sl[1].start, y+sl[0].start, peak, bck, fwhm
        
    return res0, res1  
    
def compute_fwhm_global(data, center, box):
    sl = img_box(center, data.shape, box)
    raster = data[sl]
    
    background = raster.min()
    braster = raster - background
    
    newc = center[0] - sl[0].start, center[1] - sl[1].start
    try:
        res = compute_fwhm(braster, newc)
        return res[1]+sl[1].start, res[0]+sl[0].start, res[2], res[3], res[4]
    except ValueError as error:
        _logger.warning("%s", error)
        return center[1], center[0], -99.0, -99.0, -99.0
    except Exception as error:
        _logger.warning("%s", error)
        return center[1], center[0], -199.0, -199.0, -199.0

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
        print 'compute background in annulus r1=', self.r1, 'r2=',self.r2
        bck, _ = comp_back_with_annulus(a, x, y, self.r1, self.r2)
        return bck

def enclosed_fwhm(data, x0, y0, bckestim):
    
    bck = bckestim(data, x0, y0)
    
    peak, fwhm, msg = compute_fwhm_enclosed(data, x0, y0, bck, extra=False)
    
    return x0, y0, peak, bck, fwhm, msg

def _fwhm_side_lineal(uu, vv):
    '''Compute r12 using linear interpolation.'''
    res1, = numpy.nonzero(vv < 0)
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
    
# returns x,y
def gauss_model(data, center_r):
    sl = img_box(center_r, data.shape, box=(4,4))
    raster = data[sl]
    
    # background
    background = raster.min()
    
    b_raster = raster - background
    
    new_c = center_r[0] - sl[0].start, center_r[1] - sl[1].start

    yi, xi = numpy.indices(b_raster.shape)

    g = models.Gaussian2D(amplitude=b_raster.max(), x_mean=new_c[1], y_mean=new_c[0], x_stddev=1.0, y_stddev=1.0)
    f1 = fitting.NonLinearLSQFitter()
    t = f1(g, xi, yi, b_raster)
    
    mm = t.x_mean.value + sl[1].start, t.y_mean.value+ sl[0].start, t.amplitude.value, t.x_stddev.value, t.y_stddev.value
    return mm

def pinhole_char(data, ncenters, box=4, recenter=True, maxdist=10.0):

    ibox = (box, box)
    
    # convert FITS x, y coordinates (pixel center 1)
    # to python/scipy/astropy (pixel center 0 and x,y -> y,x)
    
    npinholes = ncenters.shape[0]
    centers_py = numpy.fliplr(ncenters[:,0:2]) - 1
    
    # recentered values
    centers_r = numpy.empty_like(centers_py)
    
    # Ignore certain pinholes
    compute_mask = numpy.ones((npinholes,), dtype='bool')
    
    _logger.info('recenter pinhole coordinates')
    for idx, c in enumerate(centers_py):
        # A failsafe
        _logger.info('for pinhole %i', idx)
        _logger.info('center is x=%7.2f y=%7.2f', c[1], c[0])
        if ((c[1] > data.shape[1] - 5) or (c[1] < 5) or
            (c[0] > data.shape[0] - 5) or (c[0] < 5)):
            _logger.info('pinhole too near to the border')
            compute_mask[idx] = False
        else:
            if recenter and maxdist > 0.0:                        
                x, y, _back, msg = centering_centroid_xy(data, c[1], c[0], box=ibox, maxdist=maxdist)
                _logger.info('new center is x=%7.2f y=%7.2f', x, y)
                # Log in X,Y format                        
                _logger.debug('recenter message: %s', msg)
                centers_r[idx] = (y,x)
            else:
                centers_r[idx] = c
    
    mm0 = numpy.empty((centers_r.shape[0], 10))    
    
    # compute the FWHM without fitting
    _logger.info('compute model-free FWHM')
    fmt = 'x=%7.2f y=%7.2f peak=%6.3f fwhm_x=%6.3f fwhm_y=%6.3f'
    for idx, center in enumerate(centers_r):
        _logger.info('For pinhole %i', idx)
        if compute_mask[idx]:
            res = compute_fwhm_global(data, center, box=ibox)
            _logger.info(fmt, *res)
            mm0[idx,0:5] = res
        else:
            _logger.info('skipping')
            res = center[1], center[0], -99.0, -99.0, -99.0
        mm0[idx,0:5] = res
        
    _logger.info('compute Gaussian fitting')
    # compute the FWHM fitting a Gaussian
    for idx, center in enumerate(centers_r):
        _logger.info('For pinhole %i', idx)
        if compute_mask[idx]:
            res = gauss_model(data, center)
            fmt = 'x=%7.2f y=%7.2f peak=%6.3f stdev_x=%6.3f stdev_y=%6.3f'
            _logger.info(fmt, *res)
        else:
            _logger.info('skipping')
            res = center[1], center[0], -99.0, -99.0, -99.0
        mm0[idx, 5:8] = res[2:] 
    # Photometry in coordinates
    # x=centers_r[:,1]
    # y=centers_r[:,0]
    # with radius 2.0 pixels and 4.0 pixels
    
    apertures = [2.0, 4.0]
    _logger.info('compute photometry with aperture radii %s', apertures)
    # FIXME: aperture_circular returns values along rows, we transpose it
    mm0[:, 8:10] = photutils.aperture_circular(data, centers_r[:,1], centers_r[:,0], apertures).T
    _logger.info('done')
    # Convert coordinates to FITS
    mm0[:,0:2] += 1
    
    mm1 = numpy.empty((centers_r.shape[0], 11))
    
        # compute the FWHM without fitting
    _logger.info('compute model-free FWHM with linear fitting')
    _logger.info('compute FWHM from enclosed flux')
    fmt = 'x=%7.2f y=%7.2f peak=%6.3f fwhm_x=%6.3f fwhm_y=%6.3f'
    for idx, center in enumerate(centers_r):
        _logger.info('For pinhole %i', idx)
        if compute_mask[idx]:
            res0, res1 = compute_fwhm_global_no(data, center, box=ibox)
            #_logger.info(fmt, *res)
            mm1[idx,0:6] = res0
            mm1[idx,6:] = res1
        else:
            _logger.info('skipping')
            mm1[idx,0:2] = center[1], center[0]
            mm1[idx,2:] = -99
            
    return mm0, mm1
    
class TestPinholeRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')
    master_sky = DataProductRequirement(MasterIntensityFlat, 'Master Sky calibration')
    pinhole_nominal_positions = Requirement(CoordinateList2DType, 'Nominal positions of the pinholes')
    shift_coordinates = Parameter(True, 'Use header information to shift the pinhole positions from (0,0) to X_DTU, Y_DTU')
    box_half_size = Parameter(4, 'Half of the computation box size in pixels')
    recenter = Parameter(True, 'Recenter the pinhole coordinates')
    max_recenter_radius = Parameter(2.0, 'Maximum distance for recentering')
    

class TestPinholeRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)
    positions = Product(ArrayType)
    positions_alt = Product(ArrayType)
    DTU = Product(ArrayType)
    filter = Product(str)
    readmode = Product(str)
    IPA = Product(float)

@define_requirements(TestPinholeRecipeRequirements)
@define_result(TestPinholeRecipeResult)
class TestPinholeRecipe(BaseRecipe):

    def __init__(self):
        super(TestPinholeRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting pinhole processing')

        iinfo = gather_info_frames(rinput.obresult.frames)
        
        if iinfo:
            mode = iinfo[0]['readmode']
            if mode.lower() in EMIR_BIAS_MODES:
                use_bias = True
                _logger.info('readmode is %s, bias required', mode)
                
            else:
                use_bias = False
                _logger.info('readmode is %s, no bias required', mode)
                
        
        dark_info = gather_info_dframe(rinput.master_dark)
        flat_info = gather_info_dframe(rinput.master_flat)
        sky_info = gather_info_dframe(rinput.master_sky)

        print('images info:', iinfo)
        if use_bias:
            bias_info = gather_info_dframe(rinput.master_bias)
            print('bias info:', bias_info)
        print('dark info:', dark_info)
        print('flat info:', flat_info)
        print('sky info:', sky_info)

        # Loading calibrations
        if use_bias:
            with rinput.master_bias.open() as hdul:
                _logger.info('loading bias')
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbias)
        else:
            _logger.info('ignoring bias')
            bias_corrector = IdNode()
            
        with rinput.master_dark.open() as mdark_hdul:
            _logger.info('loading dark')
            mdark = mdark_hdul[0].data
            dark_corrector = DarkCorrector(mdark)

        exposure_corrector = DivideByExposure()

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
                fname = hdulist.filename()
                if fname:
                    _logger.info('input is %s', fname)
                else:
                    _logger.info('input is %s', hdulist)
                
                final = flow(hdulist)
                _logger.debug('output is input: %s', final is hdulist)

                
                cdata.append(final)
                
                # Files to be closed at the end
                odata.append(hdulist)
                if final is not hdulist:
                    odata.append(final)
                                 
            _logger.info("stacking %d images using 'mean'", len(cdata))
            data = combine.mean([d[0].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0][0].header.copy())

        finally:
            _logger.debug('closing images')
            for hdulist in odata:
                hdulist.close()
    
        _logger.debug('update result header')
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        
        _logger.debug('finding pinholes')
        
        try:
            filtername = hdr['FILTER']
            readmode = hdr['READMODE']
            ipa = hdr['IPA']
            xdtu = hdr['XDTU']
            ydtu = hdr['YDTU']
            zdtu = hdr['ZDTU']
        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)
        
        
        if rinput.shift_coordinates:
            #get things from header
            _logger.info('getting DTU position from header')
            _logger.info('XDTU=%6.2f YDTU=%6.2f ZDTU=%6.2f', xdtu, ydtu, zdtu)
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
        
        _logger.info('pinhole characterization')
        positions, positions_alt = pinhole_char(hdu.data, ncenters, box=rinput.box_half_size, 
                                 recenter=rinput.recenter,
                                 maxdist=rinput.max_recenter_radius)
        
        hdulist = fits.HDUList([hdu])
        
        assert hdulist[0].header['BUNIT'].lower() == 'adu/s'
        
        result = TestPinholeRecipeResult(frame=hdulist, positions=positions, 
                    positions_alt=positions_alt,
                    filter=filtername, DTU=[xdtu, ydtu, zdtu], readmode=readmode, IPA=ipa)
        return result
        
        
