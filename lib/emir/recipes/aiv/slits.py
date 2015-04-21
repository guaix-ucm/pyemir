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

import logging

import numpy
from scipy import ndimage
from scipy.ndimage.filters import median_filter
from skimage.filter import canny
import matplotlib.pyplot as plt
import matplotlib.patches

from numina.array.fwhm import compute_fwhm_2d_simple
from numina.array.utils import expand_region
from numina.core import Product, Parameter
from numina.core.requirements import ObservationResultRequirement
#
from emir.core import EmirRecipe
from emir.dataproducts import DataFrameType
from emir.dataproducts import ArrayType
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
from emir.requirements import MasterSkyRequirement

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


def get_dtur_from_header(hdr):

    # get DTU things from header
    _logger.info('getting DTU position from header')
    xdtu = hdr['XDTU']
    ydtu = hdr['YDTU']
    zdtu = hdr['ZDTU']

    # Defined even if not in the header
    xdtuf = hdr.get('XDTU_F', 1.0)
    ydtuf = hdr.get('YDTU_F', 1.0)
    xdtu0 = hdr.get('XDTU_0', 0.0)
    ydtu0 = hdr.get('YDTU_0', 0.0)
    _logger.info('XDTU=%6.2f YDTU=%6.2f ZDTU=%6.2f', xdtu, ydtu, zdtu)
    _logger.info('XDTU_F=%6.2f YDTU_F=%6.2f', xdtuf, ydtuf)
    _logger.info('XDTU_0=%6.2f YDTU_0=%6.2f', xdtu0, ydtu0)

    xdtur = (xdtu / xdtuf - xdtu0)
    ydtur = (ydtu / ydtuf - ydtu0)
    _logger.info('XDTU_R=%6.2f YDTU_R=%6.2f', xdtur, ydtur)
    dtur = [xdtur, ydtur, zdtu]
    return dtur


class TestSlitDetectionRecipe(EmirRecipe):

    # Recipe Requirements
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    median_filter_size = Parameter(4, 'Size of the median box')
    canny_sigma = Parameter(3.0, 'Sigma for the canny algorithm')
    obj_min_size = Parameter(200, 'Minimum size of the slit')
    obj_max_size = Parameter(3000, 'Maximum size of the slit')
    slit_size_ratio = Parameter(4.0, 'Minimum ratio between height and width for slits')

    # Recipe Results
    frame = Product(DataFrameType)
    slitstable = Product(ArrayType)
    DTU = Product(ArrayType)
    IPA = Product(float)

    def run(self, rinput):
        _logger.info('starting slit processing')

        _logger.info('basic image reduction')

        flow = init_filters_bdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        try:
            ipa = hdr['IPA']

            dtur = get_dtur_from_header(hdr)

        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

        _logger.debug('finding slits')

        # First, prefilter with median
        median_filter_size = rinput.median_filter_size
        canny_sigma = rinput.canny_sigma
        obj_min_size = rinput.obj_min_size
        obj_max_size = rinput.obj_max_size

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
        fill_slits = ndimage.binary_fill_holes(edges)

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
        # plt.imshow(fill_slits_clean)

        # and relabel
        _logger.debug('Label filtered objects')
        relabel_objects, nb_labels = ndimage.label(fill_slits_clean)
        _logger.debug('%d objects found after filtering', nb_labels)
        ids = range(1, nb_labels + 1)

        _logger.debug('Find regions and centers')
        regions = ndimage.find_objects(relabel_objects)
        centers = ndimage.center_of_mass(data2, labels=relabel_objects,
                                         index=ids
                                         )

        table = char_slit(data2, regions, centers,
                          slit_size_ratio=rinput.slit_size_ratio
                          )

        result = self.create_result(frame=hdulist, slitstable=table,
                                    DTU=dtur,
                                    IPA=ipa)

        return result


def char_slit(data, regions, centers, box_increase=3, slit_size_ratio=4.0):

    result = []

    for r in regions:
        _logger.debug('initial region %s', r)
        oshape = shape_of_slices(r)

        ratio = oshape[0] / oshape[1]
        if ratio < slit_size_ratio:
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

        _logger.debug('x=%f y=%f', c[1] + ref[1], c[0] + ref[0])
        _logger.debug('fwhm_x %f fwhm_y %f', fwhm_x, fwhm_y)

        colrow = ref[1] + cc + 1, ref[0] + fc + 1

        result.append([c[1] + ref[1] + 1, c[0] + ref[0] + 1, fwhm_x, fwhm_y])

        _logger.debug('Save figures slit-%d-%d', *colrow)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(datas)
        circle1 = matplotlib.patches.Circle(c[::-1], 0.6, color='r', fill=False)
        ax.add_artist(circle1)
        fig.savefig('slit-%d-%d-2d.png' % colrow)
        plt.close()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('left-rigth')
        ax.plot(datas[fc, :], 'r*-', label='%s' % colrow[0])
        ax.legend()
        fig.savefig('slit-%d-%d-lr.png' % colrow)
        plt.close()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title('top-bottom')
        ax.plot(datas[:, cc], 'r*-', label='%s' % colrow[1])
        ax.legend()
        fig.savefig('slit-%d-%d-tb.png' % colrow)
        plt.close()

    return result
