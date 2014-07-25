#
# Copyright 2008-2014 Universidad Complutense de Madrid
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

import numpy
from astropy import wcs

from numina.core import DataFrame, ObservationResult
from numina.core.recipeinout import RecipeResultAutoQC as RecipeResult

def gather_info_dframe(dataframe):
    with dataframe.open() as hdulist:
        info = gather_info_hdu(hdulist)
    return info

_meta = {
         'readmode': ('READMODE', 'undefined'),
         'bunit': ('BUNIT', 'ADU'),
         'texp': ('EXPTIME', None),
         'filter': ('FILTER', 'undefined'),
         'imagetype': ('IMGTYP', 'undefined')
         }

def gather_info_hdu(hdulist):
    

    # READMODE is STRING
    meta = {}
    meta['n_ext'] = len(hdulist)
    for key, val in _meta.items():
        meta[key] = hdulist[0].header.get(val[0], val[1])
    
    adu_s = False
    if meta['bunit'].lower() == 'adu/s':
        adu_s = True
    meta['adu_s'] = adu_s

    return meta
    
def gather_info_frames(framelist):
    iinfo = []
    for frame in framelist:
        with frame.open() as hdulist:
            iinfo.append(gather_info_hdu(hdulist))
    return iinfo

def gather_info(recipeinput):
    klass = recipeinput.__class__
    metadata = {}
    for key in klass:
        val = getattr(recipeinput, key)
        if isinstance(val, DataFrame):
            metadata[key] =  gather_info_dframe(val)
        elif isinstance(val, ObservationResult):
            metas = []
            for f in  val.frames:
                metas.append(gather_info_dframe(f))
            metadata[key] = metas
        else:
            pass
    return metadata

EMIR_BIAS_MODES = ['simple', 'bias', 'single']

def offsets_from_wcs(frames, pixref):
    '''Compute offsets between frames using WCS information.
    
    :parameter frames: sequence of FITS filenames or file descriptors
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

    with frames[0].open() as hdulist:
        wcsh = wcs.WCS(hdulist[0].header)
        skyref = wcsh.wcs_pix2sky(pixref, 1)

    for idx, frame in enumerate(frames[1:]):
        with frame.open() as hdulist:
            wcsh = wcs.WCS(hdulist[0].header)
            pixval = wcsh.wcs_sky2pix(skyref, 1)
            result[idx + 1] = -(pixval[0] - pixref[0])

    return result
