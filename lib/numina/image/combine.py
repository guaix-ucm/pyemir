#
# Copyright 2008-2009 Sergio Pascual
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

# $Id$

from __future__ import division, with_statement
import itertools
from itertools import izip
import functools

import scipy
import numpy

from numina.exceptions import Error
from numina.image._combine import internal_combine as c_internal_combine
from numina.image._combine import internal_combine2 as c_internal_combine2

__version__ = "$Revision$"


def blockgen1d(block, size):
    '''Compute 1d block intervals to be used by combine.
    
    blockgen1d computes the slices by recursively halving the initial
    interval (0, size) by 2 until its size is lesser or equal than block
    
    :param block: an integer maximum block size
    :param size: original size of the interval, it corresponds to a 0:size slice
    :return: a list of slices
    
    Example:
    
        >>> blockgen1d(512, 1024)
        [slice(0, 512, None), slice(512, 1024, None)]
        
        >>> blockgen1d(500, 1024)
        [slice(0, 256, None), slice(256, 512, None), slice(512, 768, None), slice(768, 1024, None)]
    
    '''
    def numblock(block, x):
        '''Compute recursively the numeric intervals'''
        a = x[0]
        b = x[1]
        if b - a <= block:
            return [x]
        else:
            result = []
            d = (b - a) // 2
            temp = map(numblock, [block, block], [(a, a + d), (a + d, b)])
            for i in temp:
                result.extend(i)
            return result
        
    return [slice(*l) for l in numblock(block, (0, size))]


def blockgen(blocks, shape):
    '''Generate a list of slice tuples to be used by combine.
    
    The tuples represent regions in an N-dimensional image.
    
    :param blocks: a tuple of block sizes
    :param shape: the shape of the n-dimensional array
    :return: an iterator to the list of tuples of slices
    
    Example:
        
        >>> blocks = (500, 512)
        >>> shape = (1040, 1024)
        >>> for i in blockgen(blocks, shape):
        ...     print i
        (slice(0, 260, None), slice(0, 512, None))
        (slice(0, 260, None), slice(512, 1024, None))
        (slice(260, 520, None), slice(0, 512, None))
        (slice(260, 520, None), slice(512, 1024, None))
        (slice(520, 780, None), slice(0, 512, None))
        (slice(520, 780, None), slice(512, 1024, None))
        (slice(780, 1040, None), slice(0, 512, None))
        (slice(780, 1040, None), slice(512, 1024, None))
        
    
    '''
    iterables = [blockgen1d(l, s) for (l, s) in zip(blocks, shape)]
    return itertools.product(*iterables)

def _py_internal_combine(method, nimages, nmasks, result, variance, numbers):
    for index, val in numpy.ndenumerate(result):
        values = [im[index] for im, mas in zip(nimages, nmasks) if not mas[index]]
        result[index], variance[index], numbers[index] = method(values)
        
def _py_internal_combine2(method, nimages, nmasks, noffsets,
             axis=None, dtype=None, out=None):
    return out
   
def combine_shape(shapes, offsets):
    # Computing final image size and new offsets
    sharr = numpy.asarray(shapes)
    
    offarr = numpy.asarray(offsets)        
    ucorners = offarr + sharr
    ref = offarr.min(axis=0)        
    finalshape = ucorners.max(axis=0) - ref 
    offsetsp = offarr - ref
    return (finalshape, offsetsp)
                


def combine(method, images, masks=None, offsets=None,
            result=None, variance=None, numbers=None,
            blocksize=(512, 512)):
    '''Combine images using the given combine method, with masks and offsets.
    
    Inputs and masks are a list of array objects.
    
    :param method: a string with the name of the method
    :param images: a list of 2D arrays
    :param masks: a list of 2D boolean arrays, True values are masked
    :param offsets: a list of 2-tuples
    :param result: a image where the result is stored if is not None
    :param variance: a image where the variance is stored if is not None
    :param number: a image where the number is stored if is not None
    :param blocksize: shape of the block used to combine images 
    :return: a 3-tuple with the result, the variance and the number
    :raise TypeError: if method is not callable
    
    
    Example:
    
       >>> import numpy
       >>> image = numpy.array([[1.,3.],[1., -1.4]])
       >>> inputs = [image, image + 1]
       >>> combine("mean", inputs)  #doctest: +NORMALIZE_WHITESPACE
       (array([[ 1.5,  3.5],
           [ 1.5, -0.9]]), array([[ 0.25,  0.25],
           [ 0.25,  0.25]]), array([[2, 2],
           [2, 2]], dtype=int32))
       
    '''
    
    # method should be a string
    if not isinstance(method, basestring):
        raise TypeError('method is not a string')
    
    # Check inputs
    if not images:
        raise Error("len(inputs) == 0")
        
    
    aimages = map(numpy.asanyarray, images)
    # All images have the same shape
    baseshape = aimages[0].shape
    if any(i.shape != baseshape for i in aimages[1:]):
        raise Error("Images don't have the same shape")
    
    # Additional checks
    for i in aimages:
        t = i.dtype
        if not (numpy.issubdtype(t, float) or numpy.issubdtype(t, int)):
            raise TypeError
    
    finalshape = baseshape
    offsetsp = numpy.zeros((len(images), 2))
    resize_images = True
    # Offsets
    if offsets:
        if len(images) != len(offsets):
            raise Error("len(inputs) != len(offsets)")
        
        resize_images = True
        
        finalshape, offsetsp = combine_shape(shape, offsets)
        
    amasks = [False] * len(images)
    generate_masks = False
    if masks:
        if len(images) != len(masks):
            raise Error("len(inputs) != len(masks)")
        amasks = map(functools.partial(numpy.asanyarray, dtype=numpy.bool), masks)
        for i in amasks:
            if not numpy.issubdtype(i.dtype, bool):
                raise TypeError
            if len(i.shape) != 2:
                raise TypeError
        # Error if mask and image have different shape
        if any(im.shape != ma.shape for (im, ma) in izip(aimages, amasks)):
            raise TypeError
    else:
        generate_masks = True


    nimages = aimages
    nmasks = amasks
    
    if resize_images:
        nimages = []
        nmasks = []
        for (i, o, m) in zip(aimages, offsetsp, amasks):
            newimage1 = numpy.empty(finalshape)
            newimage2 = numpy.ones(finalshape, dtype=numpy.bool) 
            pos = (slice(o[0], o[0] + i.shape[0]), slice(o[1], o[1] + i.shape[1]))
            newimage1[pos] = i
            newimage2[pos] = m
            nimages.append(newimage1)
            nmasks.append(newimage2)
    elif generate_masks:
        nmasks = []
        for (o, m) in zip(offsetsp, amasks):
            newimage2 = numpy.ones(finalshape, dtype=numpy.bool) 
            pos = (slice(o[0], o[0] + i.shape[0]), slice(o[1], o[1] + i.shape[1]))
            newimage2[pos] = m
            nmasks.append(newimage2)
        
    # Initialize results        
    if result is None:
        result = numpy.zeros(finalshape, dtype="float64")
    else:
        if result.shape != finalshape:
            raise TypeError("result has wrong shape")
    
    if variance is None:
        variance = numpy.zeros(finalshape, dtype="float64")
    else:
        if variance.shape != finalshape:
            raise TypeError("variance has wrong shape")
                
    if numbers is None:
        numbers = numpy.zeros(finalshape, dtype="int32")
    else:
        if numbers.shape != finalshape:
            raise TypeError("numbers has wrong shape")
    
    for i in blockgen(blocksize, finalshape):
        # views of the images and masks
        vnimages = [j[i] for j in nimages]
        vnmasks = [j[i] for j in nmasks]
        c_internal_combine(method, vnimages, vnmasks,
                          result[i], variance[i], numbers[i])
        
    return (result, variance, numbers)

def _combine(method, images, masks=None, offsets=None,
             dtype=None, out=None):
        # method should be a string
    if not isinstance(method, basestring):
        raise TypeError('method is not a string')
    
    # Check inputs
    if not images:
        raise Error("len(inputs) == 0")
        
    # numpy.asanyarray
    number_of_images = len(images)
    images = map(numpy.asanyarray, images)
    
    # All images have the same shape
    allshapes = [i.shape for i in images]
    baseshape = images[0].shape
    if any(shape != baseshape for shape in allshapes[1:]):
        raise Error("Images don't have the same shape")
    
    # Offsets
    if offsets:
        if len(images) != len(offsets):
            raise Error("len(inputs) != len(offsets)")
        
        finalshape, offsets = combine_shape(allshapes, offsets)
    else:
        offsets = numpy.zeros((number_of_images, 2))
        finalshape = baseshape
        
    if masks:
        if len(images) != len(masks):
            raise Error("len(inputs) != len(masks)")
    
        # Error if mask and image have different shape
        if any(imshape != ma.shape for (imshape, ma) in izip(allshapes, masks)):
            raise TypeError
    else:
        masks = [numpy.zeros(baseshape, dtype='uint8')] * number_of_images
        
    # Creating out if needed
    # We need thre numbers
    outshape = (3,) + finalshape
    if out is None:
        out = numpy.zeros(outshape, dtype='float64')
    else:
        if out.shape != outshape:
            raise TypeError("result has wrong shape")  
        
    c_internal_combine2(method, images, masks, offsets, out1=out[0], out2=out[1], out3=out[2])
    return out.astype(dtype)

def mean(data, masks=None, offsets=None,
         dtype=None, out=None):
    '''Compute the arithmetic mean along the specified axis.
    
    '''
    return _combine('mean', data, masks, offsets, dtype, out)


    
if __name__ == "__main__":
    # from numina.decorators import print_timing
  
    # Inputs
    input1 = scipy.ones((2000, 2000), dtype='int16')
    
    minputs = [input1] * 50
    minputs += [input1 * 2] * 50
    moffsets = [(1, 1), (1, 0), (0, 0), (0, 1), (-1, -1)]
    moffsets += [(1, 1), (1, 0), (0, 0), (0, 1), (-1, -1)]
    
    out = mean(minputs)
    print out[0,0,0]

    
    
