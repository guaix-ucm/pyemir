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
import scipy
from astropy.io import fits
from scipy.stats import linregress
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree as KDTree
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
# Channels are rotated!
#from emirextras.channels import CHANNELS_REG_ROT
import scipy.interpolate as itpl
import scipy.optimize as opz


from numina.core import RecipeError
from numina.core import BaseRecipe, RecipeRequirements, DataFrame
from numina.core import Requirement, Product, DataProductRequirement
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
from emir.dataproducts import CentroidsTableType

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"
            
GAUSS_FWHM_FACTOR = 2.354800

class TestPinholeRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')
    master_sky = DataProductRequirement(MasterIntensityFlat, 'Master Sky calibration')
    pinhole_nominal_positions = Requirement(CoordinateList2DType, 'Nominal positions of the pinholes')

class TestPinholeRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)
    centroids = Product(CentroidsTableType)

@define_requirements(TestPinholeRecipeRequirements)
@define_result(TestPinholeRecipeResult)
class TestPinholeRecipe(BaseRecipe):

    def __init__(self):
        super(TestPinholeRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple sky reduction')

        pepe = rinput.pinhole_nominal_positions
        print(type(pepe), pepe.dtype, pepe.shape)

        # Loading calibrations
        with rinput.master_bias.open() as hdul:
            _logger.info('loading bias')
            mbias = hdul[0].data
            bias_corrector = BiasCorrector(mbias)
            #bias_corrector = IdNode()
            
        exposure_corrector = DivideByExposure(factor=1e-3)

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

        cdata = []
        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                final = flow(hdulist)
                cdata.append(final)

            _logger.info('stacking %d images using median', len(cdata))
            data = median([d['primary'].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)

        finally:
            _logger.debug('closing images')
            for hdulist in cdata:
                hdulist.close()
            
        _logger.debug('update result header')
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdulist = fits.HDUList([hdu])


        _logger.debug('finding pinholes')
        data = hdu.data

        # To convert to scipy coordinates we have to subtract 1 and flip x<->y
        ref_centers = numpy.fliplr(rinput.pinhole_nominal_positions) - 1.0
        nl_pinholes = len(ref_centers)

        # 
        centroids = numpy.zeros((nl_pinholes, 6))

        # Compute and remove a background level
        m, sigma = background_and_sigma(data)
        _logger.info('Mean background is %7.0f', m)

        data -= m

        n_size = 5
        th = 1000
        # Search box
        boxsize = 15
        boxsize_full = 2 * boxsize

        result = find_pinholes_center_via_max(ref_centers, data, boxsize, n_size=n_size)

        for val in result:
            print val


        for idx, n_c in enumerate(ref_centers):
            n_c_i = n_c.astype('int')
            sl = (slice(n_c_i[0] - boxsize, n_c_i[0] + boxsize, None)), (slice(n_c_i[1] - boxsize, n_c_i[1] + boxsize, None))

            region = data[sl]

            data_max = filters.maximum_filter(region, n_size)
            maxima = (region == data_max)
            #data_min = filters.minimum_filter(region, neighborhood_size)
            #diff = ((data_max - data_min) > th)
            #maxima[diff == 0] = 0
            #fits.writeto('mask.fits', maxima * 1, clobber=True)

            sombra = region * maxima
            maxindex = numpy.unravel_index(sombra.argmax(), sombra.shape)
            maxvalue = sombra[maxindex]

            # Fittinh Gaussian model
            p_init = models.Gaussian2D(amplitude=maxvalue, 
             x_mean=maxindex[1], y_mean=maxindex[0], x_stddev=1.0,
             y_stddev=1.0, theta=1.0)
            f = fitting.NonLinearLSQFitter()
            x, y = numpy.mgrid[:boxsize_full, :boxsize_full]
            p = f(p_init, x, y, region)
            print 'fit center', p.x_mean, p.y_mean
            print 'fit sigma', p.x_stddev, p.y_stddev

            abs_centroid = n_c_i + [p.x_mean[0], p.y_mean[0]] - boxsize
            crow = n_c[1], n_c[0], abs_centroid[1], abs_centroid[0], p.y_stddev[0], p.x_stddev[0]
            centroids[idx, :] = crow

        # In FITS coordinates, we have to add one
        # To the first 4 columns
        centroids[:,:4] += 1

        result = find_pinholes_center_via_segmentation(ref_centers, data, boxsize, threshold=th)

        for val in result:
            print val[1]

        m = find_pinholes_via_segmentation(ref_centers, data, boxsize, threshold=th)
        print m

        result = TestPinholeRecipeResult(frame=hdulist, centroids=centroids)

        return result

def segmentation(data, threshold):
    mask = numpy.where(data > threshold, 1, 0)
    labelmask, nobj = scipy.ndimage.label(mask)
    objects = scipy.ndimage.measurements.find_objects(labelmask)
    centers = scipy.ndimage.measurements.center_of_mass(data, labels=labelmask, index=range(1, nobj+1))

    return objects, centers, labelmask

def img_box(center, shape, box):
    def slice_create(c, s, b):
        cc = int(math.floor(c + 0.5))
        l = max(0, cc - b)
        h = min(s, cc + b +1)
        return slice(l, h, None)

    return tuple(slice_create(*args) for args in zip(center, shape, box))

def find_pinholes_center_via_segmentation(ref_centers, data, boxsize, threshold):
    result = []
    for idx, n_c in enumerate(ref_centers):
        sl = img_box(n_c, data.shape, box=(boxsize,boxsize))
        region = data[sl]
        seginfo = segmentation(region, threshold)
        result.append(seginfo)
    return result

def find_pinholes_via_segmentation(ref_centers, data, boxsize, threshold):
    pinhole_fits = {}
    for idx, n_c in enumerate(ref_centers):
        sl = img_box(n_c, data.shape, box=(boxsize,boxsize))
        region = data[sl]
        seginfo = segmentation(region, threshold)
        # FIXME: choose the object with the centroinf closests to the center
        # For now, I choose the first
        n_center = seginfo[1][0]
        abs_n_center = [i + s.start for i, s in zip(n_center, sl)]
        print 'new center',n_center, abs_n_center

        # Taking a new box around the centroid
        sl = img_box(abs_n_center, data.shape, box=(15,15))
        img = data[sl]
        c = abs_n_center
        c_ff = [i - s.start for i, s in zip(c, sl)]
        try:
            params = compute_fwhm(img, c_ff)
            pinhole_fits[idx] = [params, sl, 0]
        except StandardError as error:
            # FIXME, error during fitting
            print(error)
    return pinhole_fits

def find_pinholes_center_via_max(ref_centers, data, boxsize, n_size=5):
    result = []
    for idx, n_c in enumerate(ref_centers):
        n_c_i = n_c.astype('int')
        sl = (slice(n_c_i[0] - boxsize, n_c_i[0] + boxsize, None)), (slice(n_c_i[1] - boxsize, n_c_i[1] + boxsize, None))

        region = data[sl]

        data_max = filters.maximum_filter(region, n_size)
        maxima = (region == data_max)
        #data_min = filters.minimum_filter(region, neighborhood_size)
        #diff = ((data_max - data_min) > th)
        #maxima[diff == 0] = 0
        #fits.writeto('mask.fits', maxima * 1, clobber=True)

        sombra = region * maxima
        maxindex = numpy.unravel_index(sombra.argmax(), sombra.shape)
        maxvalue = sombra[maxindex]
        result.append([maxindex, maxvalue])
    return result


def background_and_sigma(data):
    # Operations...

    # thresholding
    # Given the mask, around 91% of the image is background
    # So we can compute the value at, say 90%
    # everything below that value is background

    limit = scipy.stats.scoreatpercentile(data.ravel(), per=90)

    # Background pixels
    back = data[data < limit]
    mu = numpy.median(back)
    _logger.info('mean estimation is %f ', mu)

    back -= mu
    _logger.info('subtract mean')

    siglev = 2.0

    # Computing background Gaussian noise
    qns = 100 * scipy.stats.norm.cdf(siglev)
    pns = 100 - qns

    ls = scipy.stats.scoreatpercentile(back.flat, pns)
    hs = scipy.stats.scoreatpercentile(back.flat, qns)
    # sigma estimation
    sig = (hs - ls) / (2 * siglev)
    _logger.info('sigma estimation is %f ', sig)

    return mu, sig


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

#        plt.plot(U, V, 'r*-')
#        plt.plot([sol_l, sol_r], [0.5 * peak, 0.5 * peak])
#        plt.show()
        fwhm = sol_r - sol_l
        return fwhm

    U = X
    V = bb.ev(U, [center[1] for i in U])

    fwhm_x = compute_fwhm_1(U, V, f1, center[0])

    U = Y
    V = bb.ev([center[0] for i in U], U)

    fwhm_y = compute_fwhm_1(U, V, f2, center[1])

    return peak, center[0], center[1], fwhm_x, fwhm_y

if __name__ == '__main__':

    main()

