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

'''Recipe to detect slits in the AIV mask'''

from __future__ import division
#
import logging

import numpy

from scipy import ndimage
from scipy.ndimage.filters import median_filter

from skimage.filter import canny

from numina.array.fwhm import compute_fwhm_2d_simple
from numina.array.utils import expand_region

import matplotlib.pyplot as plt

# import math
#
# import numpy
# import scipy.interpolate as itpl
# import scipy.optimize as opz
# from astropy.modeling import models, fitting
# import photutils
#
# from numina.array.recenter import img_box, centering_centroid
#
from numina.core import BaseRecipe, RecipeRequirements, RecipeError
from numina.core import Requirement, Product, DataProductRequirement, Parameter
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement
#
from emir.core import RecipeResult
from emir.dataproducts import DataFrameType, MasterIntensityFlat
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
# from emir.dataproducts import CoordinateList2DType
# from emir.dataproducts import ArrayType

# from photutils import aperture_circular
#
# from .procedures import compute_fwhm_spline2d_fit
# from .procedures import compute_fwhm_enclosed_direct
# from .procedures import compute_fwhm_enclosed_grow
# from .procedures import compute_fwhm_simple
# from .procedures import moments
# from .procedures import AnnulusBackgroundEstimator
# from .procedures import img_box2d
from .flows import basic_processing_with_combination
from .flows import init_filters_bdfs

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"


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


class TestSlitDetectionRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = DataProductRequirement(MasterIntensityFlat,
                                        'Master Sky calibration')


class TestSlitDetectionRecipeResult(RecipeResult):
    frame = Product(DataFrameType)


@define_requirements(TestSlitDetectionRecipeRequirements)
@define_result(TestSlitDetectionRecipeResult)
class TestSlitDetectionRecipe(BaseRecipe):

    def __init__(self):
        super(TestSlitDetectionRecipe, self).__init__(
            author=_s_author,
            version="0.1.0"
            )

    def run(self, rinput):
        _logger.info('starting pinhole processing')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)
        hdr = hdulist[0].header
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')


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

        _logger.debug('finding pinholes')

        # First, prefilter with median
        median_filter_size = 4
        canny_sigma = 3
        obj_min_size = 200
        obj_max_size = 3000
        
        
        data1 = hdulist[0].data
        _logger.debug('Median filter with box %d', median_filter_size)
        data2 = median_filter(data1, size=median_filter_size)

        # Grey level image
        img_grey = normalize(data2)

        # Find edges with canny
        _logger.debug('Find edges with canny, sigma %d', canny_sigma)
        edges = canny(img_grey, sigma=canny_sigma)
        
        # Fill edges
        _logger.debug('Fill holes')
        fill_slits =  ndimage.binary_fill_holes(edges)

        _logger.debug('Label objects')
        label_objects, nb_labels = ndimage.label(fill_slits)
        _logger.debug('%d objects found', nb_labels)
        # Filter on the area of the labeled region
        # Perhaps we could ignore this filtering and
        # do it later?
        _logger.debug('Filter objects by size')
        # Sizes of regions
        sizes = numpy.bincount(label_objects.ravel())

        _logger.debug('Min size is %d', obj_min_size)
        _logger.debug('Max size is %d', obj_max_size)

        mask_sizes = (sizes > obj_min_size) & (sizes < obj_max_size)

        # Filter out regions
        nids, = numpy.where(mask_sizes)

        mm = numpy.in1d(label_objects, nids)
        mm.shape = label_objects.shape

        fill_slits_clean = numpy.where(mm, 1, 0)
        #plt.imshow(fill_slits_clean)

        # and relabel
        _logger.debug('Label filtered objects')
        relabel_objects, nb_labels = ndimage.label(fill_slits_clean)
        _logger.debug('%d objects found after filtering', nb_labels)
        ids = range(1, nb_labels + 1)

        _logger.debug('Find regions and centers')
        regions = ndimage.find_objects(relabel_objects)
        centers = ndimage.center_of_mass(data2, labels=relabel_objects, index=ids)

        char_slit(data2, regions, centers)

        result = self.create_result(frame=hdulist)

        return result

def char_slit(data, regions, centers, box_increase=3, slit_size_ratio=4.0):

    for r, c_alt in zip(regions, centers):
        print 'initial region', r
        oshape = shape_of_slices(r)

        ratio = oshape[0] / oshape[1]
        if ratio < slit_size_ratio:
            print "this is not a slit, ratio=", ratio
            continue

        print 'initial shape', oshape
        print 'ratio', ratio
        rp = expand_region(r, box_increase, box_increase,
                           start=0, stop=2048)
        print 'expanded region', rp
        ref = rp[0].start, rp[1].start
        print 'reference point', ref
    
        datas = data[rp]
        shape = datas.shape
    
        print 'data, shape', shape
        print 'orig shape', 
        print 'data, shape', shape_of_slices(rp)


        c = ndimage.center_of_mass(datas)
    
        fc = datas.shape[0] // 2
        cc = datas.shape[1] // 2
        print fc, cc, c[0], c[1]

        peak, fwhm_x, fwhm_y = compute_fwhm_2d_simple(datas, c[1], c[0])

        print 'center', 'y=',c[0] + ref[0], 'x=',c[1] +  ref[1]
        print 'center', 'y=',c_alt[0], 'x=',c_alt[1]
        print 'fwhm_x', fwhm_x    
        print 'fwhm_y', fwhm_y
    
        fig = plt.figure()
        ax = fig.add_subplot(111)  
        ax.imshow(datas)
        circle1 = plt.Circle(c[::-1], 0.6, color='r', fill=False)
        ax.add_artist(circle1)
        plt.show()
    
        plt.title('left-rigth')
        plt.plot(datas[fc,:], 'r*-', label='%s' % (ref[0] + fc + 1))
        plt.legend()
        plt.show()

        plt.title('top-bottom')
        plt.plot(datas[:,cc], 'r*-', label='%s' % (ref[1] + cc + 1))
        plt.legend()
        plt.show()
        _logger.debug('Label filtered objects')
