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

from numina.array.recenter import img_box, centering_centroid

from numina import __version__
from numina.core import BaseRecipe, RecipeRequirements, RecipeError
from numina.core import Requirement, Product, DataProductRequirement, Parameter
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement
from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import FlatFieldCorrector, SkyCorrector
from numina.flow import SerialFlow
from numina.flow.node import IdNode
from numina.array import combine

from emir.core import RecipeResult
from emir.core import EMIR_BIAS_MODES
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import DataFrameType, MasterIntensityFlat
from emir.dataproducts import CoordinateList2DType
from emir.dataproducts import ArrayType
from emir.core import gather_info
from photutils import aperture_circular

from .procedures import compute_fwhm_spline2d_fit
from .procedures import compute_fwhm_enclosed_direct
from .procedures import compute_fwhm_enclosed_grow
from .procedures import compute_fwhm_simple
from .procedures import moments
from .procedures import AnnulusBackgroundEstimator
from .procedures import img_box2d

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"
            
GAUSS_FWHM_FACTOR = 2.354800

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
    for idx, (yi,xi) in enumerate(centers_py):
        # A failsafe
        _logger.info('for pinhole %i', idx)
        _logger.info('center is x=%7.2f y=%7.2f', xi, yi)
        if ((xi > data.shape[1] - 5) or (xi < 5) or
            (yi > data.shape[0] - 5) or (yi < 5)):
            _logger.info('pinhole too near to the border')
            compute_mask[idx] = False
        else:
            if recenter and maxdist > 0.0:                        
                xc, yc, _back, msg = centering_centroid(data, xi, yi, box=ibox, maxdist=maxdist)
                _logger.info('new center is x=%7.2f y=%7.2f', xc, yc)
                # Log in X,Y format                        
                _logger.debug('recenter message: %s', msg)
                centers_r[idx] = (yc,xc)
            else:
                centers_r[idx] = (yi, xi)
    # Result 0
    mm0 = numpy.empty((centers_r.shape[0], 10))    
    mm0.fill(-99) 
    # compute the FWHM without fitting
    
    for idx, (yc, xc) in enumerate(centers_r):
        _logger.info('For pinhole %i', idx)
        mm0[idx,0:2] = xc, yc
        if compute_mask[idx]:
            _logger.info('compute model-free FWHM')
            try:
                res1 = compute_fwhm_global(data, (yc, xc), box=ibox)
                fmt1 = 'x=%7.2f y=%7.2f peak=%6.3f fwhm_x=%6.3f fwhm_y=%6.3f'
                _logger.info(fmt1, *res1)
                mm0[idx, 2:5] = res1[2:]
            except StandardError as error:
                _logger.exception("unable to obtain FWHM, %s", error)
            
            _logger.info('compute Gaussian 2Dfitting')        
            try:
                res2 = gauss_model(data, (yc, xc))
                fmt2 = 'x=%7.2f y=%7.2f peak=%6.3f stdev_x=%6.3f stdev_y=%6.3f'
                _logger.info(fmt2, *res2)
                mm0[idx, 5:8] = res2[2:]
            except StandardError as error:
                _logger.exception("unable to obtain FWHM, %s", error)
        else:
            _logger.info('skipping')
         
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
    return mm0

def pinhole_char2(data, ncenters,
        recenter=True, 
        recenter_half_box=5, 
        recenter_nloop=10,
        recenter_maxdist=10.0,
        back_buff=3,
        back_width=5,
        phot_niter=10,
        phot_rad=8):

    sigma0 = 1.0
    rad =  3 * sigma0 * GAUSS_FWHM_FACTOR
    box = recenter_half_box
    recenter_half_box = (box, box)
    
    # convert FITS x, y coordinates (pixel center 1)
    # to python/scipy/astropy (pixel center 0 and x,y -> y,x)
    
    npinholes = ncenters.shape[0]
    centers_i = ncenters[:,0:2] - 1
    
    # recentered values
    centers_r = numpy.empty_like(centers_i)

    # Ignore certain pinholes
    compute_mask = numpy.ones((npinholes,), dtype='bool')
    
    _logger.info('recenter pinhole coordinates')
    for idx, (xi, yi) in enumerate(centers_i):
        # A failsafe
        _logger.info('for pinhole %i', idx)
        _logger.info('center is x=%7.2f y=%7.2f', xi, yi)
        if ((xi > data.shape[1] - 5) or (xi < 5) or
            (yi > data.shape[0] - 5) or (yi < 5)):
            _logger.info('pinhole too near to the border')
            compute_mask[idx] = False
        else:
            if recenter and (recenter_maxdist > 0.0):                        
                xc, yc, _back, msg = centering_centroid(data, xi, yi, 
                        box=recenter_half_box, 
                        maxdist=recenter_maxdist, nloop=recenter_nloop)
                _logger.info('new center is x=%7.2f y=%7.2f', xc, yc)
                # Log in X,Y format                        
                _logger.debug('recenter message: %s', msg)
                centers_r[idx] = xc, yc
            else:
                centers_r[idx] = xi, yi
    
    # Number of results
    nresults = 32
    mm0 = numpy.empty((centers_r.shape[0], nresults))    
    mm0[:] = -99
    mm0[:,0:2] = centers_i
    mm0[:,2:4] = centers_r 

    # Fitter
    fitter = fitting.NonLinearLSQFitter()
    rplot = phot_rad
    
    for idx, (x0, y0) in enumerate(centers_r):
        _logger.info('For pinhole %i', idx)
        if not compute_mask[idx]:
            _logger.info('skipping')
            # Fill result with -99
            continue
        
        # Initial photometric radius
        rad = 3.0        
        # Loop to find better photometry radius and background annulus 
        irad = rad
        for i in range(phot_niter):
            phot_rad = rad
            # Sky background annulus
            rs1 = rad + back_buff
            rs2 = rs1 + back_width
            _logger.debug('Iter %d, annulus r1=%5.2f r2=%5.2f', i, rs1, rs2)
            bckestim = AnnulusBackgroundEstimator(r1=rs1, r2=rs2)
            
            # Crop the image to obtain the background
            sl_sky = img_box2d(x0, y0, data.shape, (rs2, rs2))
            raster_sky = data[sl_sky]
            # Logical coordinates
            xx0 = x0 - sl_sky[1].start
            yy0 = y0 - sl_sky[0].start
            # FIXME, perhaps we dont need to crop the image
            try:
                bck = bckestim(raster_sky, xx0, yy0)
                _logger.debug('Iter %d, background %f in annulus r1=%5.2f r2=%5.2f', 
                              i, bck, rs1, rs2)
            except StandardError as error:
                _logger.warning('Error in background estimation %s', error)
                break
        
            # Radius of the fit
            fit_rad = max(rplot, rad)
        
            sl = img_box2d(x0, y0, data.shape, (fit_rad, fit_rad))
            part = data[sl]
            # Logical coordinates
            xx0 = x0 - sl[1].start
            yy0 = y0 - sl[0].start
    
            Y, X = numpy.mgrid[sl]
            _logger.debug('Iter %d, radial fit', i)

            # Photometry
            D = numpy.sqrt((X-x0)**2 + (Y-y0)**2)
            phot_mask = D < fit_rad
            r1 = D[phot_mask]
            part_s = part - bck
            f1 = part_s[phot_mask]
            # Fit radial profile
            model = models.Gaussian1D(amplitude=f1.max(), mean=0, stddev=1.0)
            model.mean.fixed = True # Mean is always 0.0
            
            g1d_f = fitter(model, r1, f1, weights=(r1+1e-12)**-1)

            rpeak = g1d_f.amplitude.value
            # sometimes the fit is negative
            rsigma = abs(g1d_f.stddev.value)
            
            rfwhm = rsigma * GAUSS_FWHM_FACTOR

            rad = 2.5 * rfwhm
            _logger.debug('Iter %d, new rad is %f', i, rad)
            if abs(rad-irad) < 1e-3:
                # reached convergence
                _logger.debug('Convergence in iter %d', i)
                break
            else:
                irad = rad
        else:
            _logger.debug('no convergence in photometric radius determination')
        
        _logger.info('background %6.2f, r1 %7.2f r2 %7.2f', bck, rs1, rs2)
        mm0[idx,4:4+3] = bck, rs1, rs2
        aper_rad = rad
        flux_aper = aperture_circular(part_s, [xx0], [yy0], aper_rad)
        _logger.info('aper rad %f, aper flux %f', aper_rad, flux_aper)
        mm0[idx,7:7+2] = aper_rad, flux_aper         

        _logger.info('Radial fit, peak: %f fwhm %f', rpeak, rfwhm)
        
        try:
            dpeak, dfwhm, smsg = compute_fwhm_enclosed_direct(part_s, xx0, yy0, maxrad=fit_rad)
            _logger.info('Enclosed direct, peak: %f fwhm %f', dpeak, dfwhm)
        except StandardError as error:
            _logger.warning('Error in compute_fwhm_enclosed_direct %s', error)
            dpeak, dfwhm = -99.0, -99.0
        
        try:
            eamp, efwhm, epeak, emsg = compute_fwhm_enclosed_grow(part_s, xx0, yy0, maxrad=fit_rad)
            _logger.info('Enclosed fit, peak: %f fwhm %f', epeak, efwhm)
        except StandardError as error:
            _logger.warning('Error in compute_fwhm_enclosed_grow %s', error)
            eamp, efwhm, epeak, emsg = [-99.0] * 4 
        
        
        mm0[idx,9:9+6] = epeak, efwhm, dpeak, dfwhm, rpeak, rfwhm
    
        try:
            res_simple = compute_fwhm_simple(part_s, xx0, yy0)
            _logger.info('Simple, peak: %f fwhm x %f fwhm %f', *res_simple)
            mm0[idx,15:15+3] = res_simple
        except StandardError as error:
            _logger.warning('Error in compute_fwhm_simple %s', error)
            mm0[idx,15:15+3] = -99.0

        try:
            res_spline = compute_fwhm_spline2d_fit(part_s, xx0, yy0)
            _logger.info('Spline, peak: %f fwhm x %f fwhm %f', *res_spline)
            mm0[idx,18:18+3] = res_spline
        except StandardError as error:
            _logger.warning('Error in compute_fwhm_spline2d_fit %s', error)
            mm0[idx,18:18+3] = -99.0
        
        # Bidimensional fit
        # Fit in a smaller box
        fit2d_rad = int(math.ceil(fit_rad))
    
        fit2d_half_box = (fit2d_rad, fit2d_rad)
        sl1 = img_box2d(x0, y0, data.shape, fit2d_half_box)
        
        part1 = data[sl1]
        Y1, X1 = numpy.mgrid[sl1]
        
        g2d = models.Gaussian2D(amplitude=rpeak, x_mean=x0, y_mean=y0, 
                                x_stddev=1.0, y_stddev=1.0)
        g2d_f = fitter(g2d, X1, Y1, part1 - bck)
        
        res_gauss2d = (g2d_f.amplitude.value, 
                        g2d_f.x_mean.value + 1, # FITS coordinates 
                        g2d_f.y_mean.value + 1, # FITS coordinates
                        g2d_f.x_stddev.value * GAUSS_FWHM_FACTOR,
                        g2d_f.y_stddev.value * GAUSS_FWHM_FACTOR,
                        g2d_f.theta.value)

        _logger.info('Gauss2d, %s', res_gauss2d)
        mm0[idx,21:21+6] = res_gauss2d
        # Moments
        moments_half_box = fit2d_half_box
        res_moments = moments(data, x0, y0, moments_half_box)
        _logger.info('Mxx %f Myy %f Mxy %f e %f pa %f', *res_moments)
        
        mm0[idx,27:27+5] = res_moments

    # FITS coordinates
    mm0[:,:4] += 1

    return mm0
    
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
    frame = Product(DataFrameType)
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
        
        meta = gather_info(rinput)
        iinfo = meta['obresult']
        
        if iinfo:
            mode = iinfo[0]['readmode']
            if mode.lower() in EMIR_BIAS_MODES:
                use_bias = True
                _logger.info('readmode is %s, bias required', mode)
                
            else:
                use_bias = False
                _logger.info('readmode is %s, no bias required', mode)
                
        dark_info = meta['master_dark']
        flat_info = meta['master_flat']
        sky_info = meta['master_sky']

        print('images info:', iinfo)
        if use_bias:
            bias_info = meta['master_bias']
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


        with rinput.master_flat.open() as mflat_hdul:
            _logger.info('loading intensity flat')
            mflat = mflat_hdul[0].data
            flat_corrector = FlatFieldCorrector(mflat)

        with rinput.master_sky.open() as msky_hdul:
            _logger.info('loading sky')
            msky = msky_hdul[0].data
            sky_corrector = SkyCorrector(msky)

        flow = SerialFlow([bias_corrector,
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
        positions = pinhole_char(hdu.data, 
                 ncenters, box=rinput.box_half_size, 
                 recenter=rinput.recenter,
                 maxdist=rinput.max_recenter_radius)
        
        _logger.info('alternate pinhole characterization')
        positions_alt = pinhole_char2(hdu.data, ncenters, 
            recenter=rinput.recenter,
            recenter_half_box=rinput.box_half_size, 
            recenter_maxdist=rinput.max_recenter_radius)
        
        hdulist = fits.HDUList([hdu])
        
        result = TestPinholeRecipeResult(frame=hdulist, positions=positions, 
                    positions_alt=positions_alt,
                    filter=filtername, DTU=[xdtu, ydtu, zdtu], readmode=readmode, IPA=ipa)
        return result
        
        
