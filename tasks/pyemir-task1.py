#!/usr/bin/python

#
# Copyright 2008-2012 Universidad Complutense de Madrid
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

from __future__ import print_function
from __future__ import division

import sys
import argparse

import numpy
import pyfits

def preproc_single(data):
    return data

def preproc_cds(data):
    return data[...,1] - data[...,0]

taskname = 'Imaging flats'
taskidx = '2.1.1'
description = "Execute '{}', task code {}".format(taskname, taskidx)

parser = argparse.ArgumentParser(description=description)
parser.add_argument("flatlist", help='list of on-flats and off-flats')
parser.add_argument("--output", help='FITS output', default='mean.fits')
args = parser.parse_args()

SINGLE, CDS = 0,1

ll = []

with open(args.flatlist) as fd:
    for l in fd:
        ll.append(l[:-1])

n2 = len(ll)
hp = n2 // 2

if 2 * hp != n2:
    print('Error, number of flat images is odd')
    sys.exit(1)

# All the images are the same shape:
with pyfits.open(ll[0]) as hdul:
    shape2d = hdul[0].shape
    print('Shape of first image is {0}'.format(shape2d))
    ndims = len(shape2d)
    if ndims == 2:
        print('Images are 2D')
        print('Asuming mode is SINGLE')
        mode = SINGLE
        preproc = preproc_single
    elif ndims == 3 and shape2d[2] == 2:
        print('Images are 3D, 3-axis is 2')
        print('Asuming mode is CDS')
        mode = CDS
        preproc = preproc_cds
    else:
        print('Errror: image has dimension={} and shape={}', ndims, shape2d)
        print('I dont\' know what to do')
        sys.exit(1)
        
    shape3d = shape2d[0], shape2d[1], hp

final = numpy.empty(shape=shape3d, dtype='int32')

for idx, (conf, sinf) in enumerate(zip(ll[:hp], ll[hp:])):
    print('Subtracting: {0} - {1}'.format(conf, sinf))
    data0 =  pyfits.getdata(conf)
    data1 =  pyfits.getdata(sinf)
    final[..., idx] = preproc(data0) - preproc(data1)

print('Computing the mean of {0} layers'.format(hp))
mdata = final.mean(axis=-1, dtype='float32')

mhdu = pyfits.PrimaryHDU(mdata)
print('Writing mean to {0}'.format(args.output))
mhdu.writeto(args.output, clobber=True)

regions = [
    (slice(0,1024, None), slice(0,1024,None)),
    (slice(0,1024, None), slice(1024,2048,None)),
    (slice(1024,2048, None), slice(0,1024,None)),
    (slice(1024,2048, None), slice(1024,2048,None))]

for idx, region in enumerate(regions):
    rr = mdata[region]
    print('Quadrant {}, mean={}, rms={}'.format(idx + 1, rr.mean(), rr.std()))

