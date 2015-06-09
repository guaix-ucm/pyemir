# Version 29 May 2015
#------------------------------------------------------------------------------

from __future__ import division
from __future__ import print_function

import numpy as np

#------------------------------------------------------------------------------

def zscale(image, factor = 0.25, LDEBUG = False):
    '''Compute z1 and z2 cuts in a similar way to Iraf.
    
    If the total number of pixels is less than 10, the function simply returns
    the minimum and the maximum values.
    '''

    npixels = image.size

    if npixels < 10:
        z1 = np.min(image)
        z2 = np.max(image)
    else:
        fnpixels = float(npixels)
        q000, q375, q500, q625, q1000 = \
          np.percentile(image,[00.0, 37.5, 50.0, 62.5, 100.0])
        zslope = (q625-q375)/(0.25*fnpixels)
        z1 = q500-(zslope*fnpixels/2)/factor
        z1 = max(z1, q000)
        z2 = q500+(zslope*fnpixels/2)/factor
        z2 = min(z2, q1000)
    
    if LDEBUG:
        print('>>> z1, z2:',z1,z2)

    return z1, z2

