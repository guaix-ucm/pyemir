#
# Copyright 2008 Sergio Pascual
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

import numpy
from scipy import maximum, minimum

def subarray_match(shape, ref, sshape, sref = (0,0)):
  '''Matches two arrays, given the different shapes given by shape and sshape, and a reference point ref
  It returnes a tuple of slices, that can be passed to the two images as indexes
  
  Example: 
  im = numpy.zeros((1000,1000))
  sim = numpy.ones((40, 40))
  i,j = subarray_match(im.shape, [20, 23], sim.shape)
  im[i] = 2 * sim[j]
  
  '''
  # Reference point in im
  ref1 = numpy.array(ref, dtype='int')  
  # Ref2 is the reference point in sim  
  # The pixel [0,0] of the array by default
  ref2 = numpy.array(sref)   
  offset = ref1 - ref2    
  urc1 = minimum(offset + numpy.array(sshape) - 1, numpy.array(shape) - 1)
  blc1 = maximum(offset, 0)
  urc2 = urc1 - offset
  blc2 = blc1 - offset
  return (slice(blc1[0],urc1[0] + 1), slice(blc1[1],urc1[1] + 1)), (slice(blc2[0],urc2[0] + 1),slice(blc2[1],urc2[1] + 1))
