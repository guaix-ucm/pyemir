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

import logging

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

from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import FrameDataProduct, MasterIntensityFlat
from emir.dataproducts import DarkCurrentValue, CoordinateList2DType
from emir.dataproducts import CentroidsTableType

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"
            
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
        if rinput.master_bias:
            _logger.info('loading bias')
            with rinput.master_bias.open() as hdul:
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbias)
        else:
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

        flow = SerialFlow([bias_corrector, dark_corrector, flat_corrector, sky_corrector])

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

        centroids = numpy.zeros((nl_pinholes, 6))
        m, sigma = threshold(data)

        data -= m

        neighborhood_size = 5
        th = 2000
        boxsize = 15
        boxsize_full = 2 * boxsize

        for idx, n_c in enumerate(ref_centers):
            n_c_i = n_c.astype('int')
            sl = (slice(n_c_i[0] - boxsize, n_c_i[0] + boxsize, None)), (slice(n_c_i[1] - boxsize, n_c_i[1] + boxsize, None))

            region = data[sl]
            #fits.writeto('a.fits', region, clobber=True)

            data_max = filters.maximum_filter(region, neighborhood_size)
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

            if False:
                plt.imshow(region)
                plt.plot(maxindex[1], maxindex[0], 'ro')
                plt.show()
                plt.imshow(p(x,y))
                plt.plot(maxindex[1], maxindex[0], 'ro')
                plt.show()
                plt.imshow(region-p(x,y))
                plt.plot(maxindex[1], maxindex[0], 'ro')
                plt.show()
            #    plt.savefig('/tmp/result.png', bbox_inches = 'tight')

        #image_fits = find_pinholes(ref_centers, centers, objects, data, prefix='pinhole-large-%d.png')

        # In FITS coordinates, we have to add one
        centroids += 1

        result = TestPinholeRecipeResult(frame=hdulist, centroids=centroids)

        return result

GAUSS_FWHM_FACTOR = 2.354800

def gaussian_checks(img):
#    params, fit = fitgaussian(img)
    params = fitgaussian(img)
    fit = gaussian(*params)
    return params, fit

def find_pinholes(ref_centers, centers, objects, data, prefix='pinhole-%d.png'):
    kdtree = KDTree(ref_centers)
    d, i = kdtree.query(centers, distance_upper_bound=2)

    pinhole_fits = {}

    for idx, (dis, el) in enumerate(zip(d, i)):
        if numpy.isfinite(dis):
#            print('idx', idx, dis, 'ref-pinhole=', el, centers[idx], ref_centers[el])
#            print('object-idx=', idx, 'is paired with ref-pinhole=', el)
            sl = objects[idx]
            img = data[sl]
#            pyfits.writeto('mytestimg.fits', img, clobber=True)

#            params, fit = fitgaussian(img)
            params, fit = gaussian_checks(img)
            pinhole_fits[el] = params
#            print(params)
#            print(params[3], params[4])
#            print(GAUSS_FWHM_FACTOR * params[3], GAUSS_FWHM_FACTOR * params[4])
            plt.figure()
            plt.imshow(img)
#            m = numpy.indices(img.shape)
#            m.shape = (2, -1)
#            plt.contour(fit(m.T).reshape(img.shape), cmap=cm.copper)
            plt.contour(fit(*numpy.indices(img.shape)), cmap=cm.copper)
            plt.savefig(prefix % idx)
            plt.close()

    return pinhole_fits

def segmentation(data, _logger):
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

    threshold = 6 * sig
    _logger.info('threshold is %f ', threshold)

    mask = numpy.where(data > threshold, 1, 0)
    fits.writeto('gmask', mask, clobber=True)
    labelmask, nobj = scipy.ndimage.label(mask)

    objects = scipy.ndimage.measurements.find_objects(labelmask)
    centers = scipy.ndimage.measurements.center_of_mass(data, labels=labelmask, index=range(1, nobj+1))

    return objects, centers

def threshold(data):
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

if __name__ == '__main__':

    main()
if __name__ == '__main__':

    main()
