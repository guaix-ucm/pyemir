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

'''Preprocessing EMIR readout modes'''

import pyfits
import numpy

def fits_wrapper(frame):
    if isinstance(frame, basestring):
        return pyfits.open(frame)
    elif isinstance(frame, pyfits.HDUList):
        return frame
    else:
        raise TypeError
        
    
def preprocess(frame):
    with fits_wrapper(frame) as hdu:
        header = hdu[0].header
        
        if header.get('READPROC'):
            return frame
        
        readmode = header.get('READMODE', "") == 'SINGLE'
        
        if readmode == 'SINGLE': 
            pass
        elif readmode == 'CDS':
            pass
        elif readmode == 'FOWLER':
            pass
        elif readmode == 'RAMP':
            pass
        else:
            raise ValueError
    pass

def preprocess_single(frame, saturation=65536, badpixels=None, blank=0):
    pass
    

def preprocess_cds(frame):
    
    samp1 = frame[0].data[...,0]
    samp2 = frame[0].data[...,1]
    final = samp2 - samp1
    frame[0].data = final
    frame[0].header['READPROC'] = True
    return frame

def preprocess_fowler(frame, saturation=65536, badpixels=None, blank=0):

    if frame[0].header['readmode'] != 'FOWLER':
        raise ValueError('Frame is not in Fowler mode')


    # We are assuming here that the 3rd axis is the number of frames
    shape = frame[0].shape
    # This construction could be handled inside loopover_fowler
    if badpixels is None:
        badpixels = numpy.zeros(shape[0:-1], dtype='>u1')
        
    hsize = shape[2] // 2
    if 2 * hsize != shape[2]:
        raise ValueError('Number of samples must be an even number')

    img, var, nmap, mask = loopover_fowler(frame[0].data, badpixels, saturation, hsize, blank=blank)

    frame[0].data = img
    frame[0].header['READPROC'] = True
    varhdu = pyfits.ImageHDU(var)
    varhdu.update_ext_name('Variance')
    frame.append(varhdu)
    nmap = pyfits.ImageHDU(nmap)
    nmap.update_ext_name('MAP')
    frame.append(nmap)
    nmask = pyfits.ImageHDU(mask)
    nmask.update_ext_name('MASK')
    frame.append(nmask)
    return frame

def axis_fowler(data, badpix, img, var, nmap, mask, hsize, saturation, blank=0):
    '''Apply Fowler processing to a series of data.'''
    MASK_SATURATION = 3 
    MASK_GOOD = 0
     
    if badpix[0] != MASK_GOOD:
        img[...] = blank
        var[...] = blank
        mask[...] = badpix[0]
    else:
        mm = numpy.asarray([(b - a) for a,b in zip(data[:hsize],data[hsize:]) if b < saturation and a < saturation])
        npix = len(mm)
        nmap[...] = npix
        if npix == 0:
            img[...] = blank
            var[...] = blank
            mask[...] = MASK_SATURATION
        elif npix == 1:
            img[...] = mm.mean()
            var[...] = blank
            mask[...] = MASK_GOOD
        else:
            img[...] = mm.mean()
            var[...] = mm.var() / npix
            mask[...] = MASK_GOOD


def loopover_fowler(data, badpixels, saturation, hsize, blank=0):
    '''Loop over the 3d array applying Fowler processing.'''
    imgfin = None
    varfin = None

#    imgfin = numpy.empty(badpixels.shape, dtype='>i2') # int16, bigendian
#    varfin = numpy.empty(badpixels.shape, dtype='>i2') # int16, bigendian

    nmask = numpy.empty(badpixels.shape, dtype='>u1') # uint8, bigendian
    npixmask = numpy.empty(badpixels.shape, dtype='>u1') # uint8, bigendian

    it = numpy.nditer([data, badpixels, imgfin, varfin, npixmask, nmask], 
                flags=['reduce_ok', 'external_loop',
                    'buffered', 'delay_bufalloc'],
                    op_flags=[['readonly'], ['readonly', 'no_broadcast'], 
                            ['readwrite', 'allocate'], 
                            ['readwrite', 'allocate'],
                            ['readwrite', 'allocate'],
                            ['readwrite', 'allocate'], 
                            ],
                    op_axes=[None,
                            [0,1,-1], 
                            [0,1,-1],
                            [0,1,-1], 
                            [0,1,-1], 
                            [0,1,-1]
                           ])
    for i in range(2, 6):
        it.operands[i][...] = 0
    it.reset()

    for x, badpix, img, var, nmap, mask in it:
        axis_fowler(x, badpix, img, var, nmap, mask, hsize, saturation, blank=blank)

    # Building final frame
    return tuple(it.operands[i] for i in range(2, 6))


def preprocess_ramp(frame):
    pass
