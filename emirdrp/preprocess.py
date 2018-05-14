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

"""Preprocessing EMIR readout modes"""

from __future__ import division

import six

from astropy.io import fits

from numina.array import ramp_array, fowler_array

from .core import EMIR_READ_MODES

PREPROC_KEY = 'READPROC'
PREPROC_VAL = True


class ReadModeGuessing(object):
    def __init__(self, mode, info=None):
        self.mode = mode
        self.info = info


def image_readmode(hdulist, default=None):
    header = hdulist[0].header
    if 'READMODE' in header:
        p_readmode = header['READMODE'].lower()
        if p_readmode in EMIR_READ_MODES:
            return ReadModeGuessing(mode=p_readmode,
                                    info={'source': 'keyword'}
                                    )
    # Using heuristics
    shape = hdulist[0].shape
    if len(shape) == 2:
        # A 2D image, mode is single
        return ReadModeGuessing(mode='single', info={'source': 'heuristics'})
    else:
        nd = min(shape)
        if nd == 2:
            # A NXNX2 image
            # mide is cds
            return ReadModeGuessing(mode='cds', info={'source': 'heuristics'})
    # Insufficient information
    if default:
        return ReadModeGuessing(mode=default, info={'source': 'default'})
    else:
        return None


def preprocess_single(hdulist):
    return hdulist


def preprocess_cds(hdulist):
    # CDS is Fowler with just one pair of reads
    return preprocess_fowler(hdulist)


def preprocess_fowler(hdulist):
    hdulist[0].header.update(PREPROC_KEY, PREPROC_VAL)
    # We need:
    tint = 0.0  # Integration time (from first read to last read)
    ts = 0.0  # Time between samples
    gain = 1.0  # Detector gain (number)
    ron = 1.0  # Detector RON (number)
    # A master badpixel mask
    cube = hdulist[0].data

    res, var, npix, mask = fowler_array(cube, ti=tint, ts=ts,
                                        gain=gain, ron=ron,
                                        badpixels=None,
                                        dtype='float32',
                                        saturation=55000.0
                                        )
    hdulist[0].data = res
    varhdu = fits.ImageHDU(var)
    varhdu.update_ext_name('VARIANCE')
    hdulist.append(varhdu)
    nmap = fits.ImageHDU(npix)
    nmap.update_ext_name('MAP')
    hdulist.append(nmap)
    nmask_hdu = fits.ImageHDU(mask)
    nmask_hdu.update_ext_name('MASK')
    hdulist.append(nmask_hdu)
    return hdulist


def preprocess_ramp(hdulist):
    hdulist[0].header.update(PREPROC_KEY, PREPROC_VAL)
    cube = hdulist[0].data
    # We need
    ti = 0.0  # Integration time
    gain = 1.0
    ron = 1.0
    rslt = ramp_array(cube, ti, gain=gain,
                      ron=ron, badpixels=None,
                      dtype='float32', saturation=55000.0
                      )
    result, var, npix, mask = rslt

    hdulist[0].data = result
    varhdu = fits.ImageHDU(var)
    varhdu.update_ext_name('VARIANCE')
    hdulist.append(varhdu)

    nmap = fits.ImageHDU(npix)
    nmap.update_ext_name('MAP')
    hdulist.append(nmap)
    nmask = fits.ImageHDU(mask)
    nmask.update_ext_name('MASK')
    hdulist.append(nmask)
    return hdulist


def fits_wrapper(frame):
    if isinstance(frame, six.string_types):
        return fits.open(frame)
    elif isinstance(frame, fits.HDUList):
        return frame
    else:
        raise TypeError


def preprocess(input_, output):
    with fits_wrapper(input_) as hdulist:
        header = hdulist[0].header
        if 'PREPROC' in header:
            # if the image is preprocessed, do nothing
            if input != output:
                hdulist.writeto(output, overwrite=True)
            return
        # determine the READ mode
        guess = image_readmode(hdulist, 'single')
        if guess is None:
            # We have a problem here
            return

        if guess.mode == 'single':
            hduproc = preprocess_single(hdulist)
        elif guess.mode == 'cds':
            hduproc = preprocess_cds(hdulist)
        elif guess.mode == 'fowler':
            hduproc = preprocess_fowler(hdulist)
        elif guess.mode == 'ramp':
            pass
        else:
            hduproc = preprocess_single(hdulist)

        hduproc.writeto(output, overwrite=True)
