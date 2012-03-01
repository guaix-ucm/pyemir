#
# Copyright 2011-2012 Universidad Complutense de Madrid
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

'''

Routines shared by image mode recipes

'''

import numpy
import pyfits
import pywcs

def name_redimensioned_images(label, iteration, ext='.fits'):
    dn = '%s_r%s' % (label, ext)
    mn = '%s_mr%s' % (label, ext)
    return dn, mn

def name_object_mask(label, iteration, ext='.fits'):
    return '%s_mro_i%01d%s' % (label, iteration, ext)

def name_skybackground(label, iteration, ext='.fits'):
    dn = '%s_sky_i%01d%s' % (label, iteration, ext)
    return dn

def name_skybackgroundmask(label, iteration, ext='.fits'):
    dn = '%s_skymask_i%01d%s' % (label, iteration, ext)
    return dn

def name_skysub_proc(label, iteration, ext='.fits'):
    dn = '%s_rfs_i%01d%s' % (label, iteration, ext)
    return dn

def name_skyflat(label, iteration, ext='.fits'):
    dn = 'superflat_%s_i%01d%s' % (label, iteration, ext)
    return dn

def name_skyflat_proc(label, iteration, ext='.fits'):
    dn = '%s_rf_i%01d%s' % (label, iteration, ext)
    return dn

def name_segmask(iteration, ext='.fits'):
    return "check_i%01d%s" % (iteration, ext)


def offsets_from_wcs(frames, pixref):
    '''Compute offsets between frames using WCS information.
    
    :parameter frames: sequence of FITS filenames or file descritors
    :parameter pixref: numpy array used as reference pixel
    
    The sky world coordinates are computed on *pixref* using
    the WCS of the first frame in the sequence. Then, the
    pixel coordinates of the reference sky world-coordinates 
    are computed for the rest of the frames.
    
    The results is a numpy array with the difference between the
    computed pixel value and the reference pixel. The first line
    of the array is [0, 0], being the offset from the first image
    to itself. 
    
    '''
    
    result = numpy.zeros((len(frames), pixref.shape[1]))

    with pyfits.open(frames[0]) as hdulist:
        wcs = pywcs.WCS(hdulist[0].header)
        skyref = wcs.wcs_pix2sky(pixref, 1)

    result[0] = pixref[0] - pixref[0]

    for idx, img in enumerate(frames[1:]):
        with pyfits.open(img) as hdulist:
            wcs = pywcs.WCS(hdulist[0].header)
            pixval = wcs.wcs_sky2pix(skyref, 1)
            result[idx + 1] = pixval[0] - pixref[0]

    return result

