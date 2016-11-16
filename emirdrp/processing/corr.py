#
# Copyright 2016 Universidad Complutense de Madrid
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


import numpy
import scipy.signal
import numina.array.imsurfit as imsurfit


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


def offsets_from_crosscor(arrs, region, refine=True, refine_box=3, order='ij'):

    # allowed values for order
    if order not in ['xy', 'ij']:
        raise ValueError("'order' must be either 'ij' or 'xy'")

    result = numpy.zeros((len(arrs), 2))
    ref_array = arrs[0]
    shape = ref_array.shape
    data1 = ref_array[region]
    d1 = standarize(data1)
    dcenter = numpy.asarray(d1.shape) // 2
    # Correlate
    for idx, arr in enumerate(arrs[1:], 1):

        data2 = arr[region]
        d2 = standarize(data2)
        corr = scipy.signal.correlate2d(d1, d2, mode='same', boundary='symm')
        # normalize
        corr /= corr.max()
        # fits.writeto('corr_%d.fits' % idx, corr, clobber=True)
        # Find peak in cross-cor
        maxindex = numpy.unravel_index(corr.argmax(), corr.shape)
        baseoff = dcenter - numpy.asarray(maxindex)
        # Pixel (0,0) in reference corresponds to baseoff in image
        # print("Pixel (0,0) in reference corresponds to %s in image" % baseoff)

        if refine:
            # Refine to subpixel
            # Fit a 2D surface to the peak of the crosscorr
            region_ref = utils.image_box(maxindex, shape, box=(refine_box, refine_box))

            coeffs, = imsurfit.imsurfit(corr[region_ref], order=2)
            # coefss are a + b *x + c * y + d*x**2 + e * x* y + f * y**2

            # get peak from coeffs
            xm, ym = vertex_of_quadratic(coeffs)
            # xm ym are offsets

            if abs(xm) > 1 or abs(ym) > 1:
                # probably bad fit
                # dont apply
                final = maxindex
            else:
                final = maxindex + numpy.asarray([ym, xm])
        else:
            final = maxindex

        refoff = dcenter - final
        # Pixel (0,0) in reference corresponds to baseoff in image
        # print("Pixel (0,0) in reference corresponds to %s in image" % refoff)
        # result xy
        if order == 'xy':
            result[idx] = refoff[::-1]
        else:
            result[idx] = refoff

    return result


if __name__ == '__main__':

    import scipy.ndimage
    import numpy.random
    import numina.array.utils as utils

    def generate_img():
        result = []
        shape = (1000, 1000)
        off = [(500, 500), (480, 490), (505, 500), (500, 500)]

        for idx, o in enumerate(off):
            data = numpy.zeros(shape) + 1000
            data[o] = 50000
            data = scipy.ndimage.gaussian_filter(data, 2.0)
            data = numpy.random.normal(data, 30.0)
            result.append(data)

        return result

    arrs = generate_img()

    shape = arrs[0].shape
    xref_cross = shape[1] // 2
    yref_cross = shape[0] // 2
    box = 50
    region = utils.image_box2d(xref_cross, yref_cross, shape, (box, box))
    print(offsets_from_crosscor(arrs, region, refine=False))
