#
# Copyright 2016-2017 Universidad Complutense de Madrid
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

"""Offsets from cross-correlation"""


from __future__  import print_function

import logging

import numpy
import numpy.linalg
import scipy.signal
import numina.array.imsurfit as imsurfit
import numina.array.utils as utils
import numina.array.stats as s


_logger = logging.getLogger('numina.recipes.emir')


def standarize(arr):
    m = arr.mean()
    s = arr.std()
    return (arr - m) / s


def vertex_of_quadratic(coeffs):
    C = coeffs[1]
    D = coeffs[2]
    A = coeffs[3]
    E = coeffs[4]
    B = coeffs[5]

    det = 4 * A * B - E**2
    if det <= 0:
        raise ValueError('quadratic has no maximum')

    xm = -(2*B*C- D*E) / det
    ym = -(2*A*D- C*E) / det
    return xm, ym


def filter_region(arr, level=4):
    median = numpy.median(arr)
    std = s.robust_std(arr, debug=False)

    return numpy.where(arr >= median + level * std, arr, 0.0)


def offsets_from_crosscor(arrs, region, refine=True, refine_box=3, order='ij'):
    # import astropy.io.fits as fits
    # allowed values for order
    if order not in ['xy', 'ij']:
        raise ValueError("'order' must be either 'ij' or 'xy'")

    result = numpy.zeros((len(arrs), 2))
    ref_array = arrs[0]

    for idx, arr in enumerate(arrs[1:], 1):

        refoff = offset_from_crosscor(ref_array, arr, region,
                                      refine=refine,
                                      refine_box=refine_box,
                                      order=order
                                      )
        result[idx] = refoff

    return result


def offsets_from_crosscor_regions(arrs, regions, refine=True, refine_box=3, order='ij', tol=0.5):
    # import astropy.io.fits as fits
    # allowed values for order
    if order not in ['xy', 'ij']:
        raise ValueError("'order' must be either 'ij' or 'xy'")

    result = numpy.zeros((len(arrs), 2))
    ref_array = arrs[0]

    for idx, arr in enumerate(arrs[1:], 1):

        refoff = offset_from_crosscor_regions(
            ref_array, arr, regions,
            refine=refine,
            refine_box=refine_box,
            order=order,
            tol=tol
        )
        result[idx] = refoff

    return result


def offset_from_crosscor(arr0, arr1, region, refine=True, refine_box=3, order='ij'):
    # import astropy.io.fits as fits
    # allowed values for order
    if order not in ['xy', 'ij']:
        raise ValueError("'order' must be either 'ij' or 'xy'")

    ref_array = arr0
    shape = ref_array.shape
    data1 = filter_region(ref_array[region])
    d1 = standarize(data1)
    # fits.writeto('cutout_%d.fits' % 0, d1, overwrite=True)

    dcenter = numpy.asarray(d1.shape) // 2
    # Correlate

    data2 = filter_region(arr1[region])
    d2 = standarize(data2)
    # fits.writeto('cutout_%d.fits' % idx, d2, overwrite=True)
    #corr = scipy.signal.correlate2d(d1, d2, mode='same', boundary='fill', fillvalue=fillvalue)
    # correlation is equivalent to convolution with inverted image
    corr = scipy.signal.fftconvolve(d1, d2[::-1, ::-1], mode='same')
    # normalize
    corr /= corr.max()
    # fits.writeto('corr_%d.fits' % idx, corr, overwrite=True)
    # fits.writeto('corr2_%d.fits' % idx, corr2 / corr2.max(), overwrite=True)
    # Find peak in cross-cor
    maxindex = numpy.unravel_index(corr.argmax(), corr.shape)

    # Check the peak is above n times the background
    median_corr = numpy.median(corr)
    std_corr = numpy.std(corr)
    peakvalue = corr[maxindex]
    threshold = median_corr + 5 * std_corr
    _logger.debug('The peak value is %f, threshold %f', peakvalue, threshold)
    if peakvalue < threshold:
        raise ValueError('Peak below threshold')

    #baseoff = dcenter - numpy.asarray(maxindex)
    # Pixel (0,0) in reference corresponds to baseoff in image
    # print("Pixel (0,0) in reference corresponds to %s in image" % baseoff)

    if refine:
        # Refine to subpixel
        # Fit a 2D surface to the peak of the crosscorr
        region_ref = utils.image_box(maxindex, shape, box=(refine_box, refine_box))

        coeffs, = imsurfit.imsurfit(corr[region_ref], order=2)
        # coefss are a + b *x + c * y + d*x**2 + e * x* y + f * y**2

        try:
            # get peak from coeffs
            xm, ym = vertex_of_quadratic(coeffs)
            # xm ym are offsets

            if abs(xm) > 1 or abs(ym) > 1:
                # probably bad fit
                # dont apply
                final = maxindex
            else:
                final = maxindex + numpy.asarray([ym, xm])
        except ValueError as error:
            _logger.debug('Error fitting peak, %s', error)
            final = maxindex
    else:
        final = maxindex

    refoff = dcenter - final
    # Pixel (0,0) in reference corresponds to baseoff in image
    # print("Pixel (0,0) in reference corresponds to %s in image" % refoff)
    # result xy
    if order == 'xy':
        return refoff[::-1]
    else:
        return refoff


def offset_from_crosscor_regions(arr0, arr1, regions, refine=True, refine_box=3, order='ij', tol=0.5):

    values = []
    for region in regions:
        try:
            res = offset_from_crosscor(arr0, arr1, region,
                                       refine=refine,
                                       refine_box=refine_box,
                                       order=order)
            values.append(res)
        except ValueError as error:
            _logger.debug('error in offset_from_crosscor_regions: %s', error)

    if len(values) == 0:
        raise ValueError('No measurements to compute offset in any region')
    values = numpy.array(values)
    _logger.debug('offsets from cross-correlation:\n%s', values)
    med_c = numpy.median(values, axis=0)

    dis = numpy.linalg.norm(values - med_c, axis=1)
    # filter values with dis > tol
    mask = dis <= tol
    # Compute mean
    if numpy.any(mask):
        return values[mask].mean(axis=0)
    else:
        raise ValueError('No measurements within tolerance in any region')
