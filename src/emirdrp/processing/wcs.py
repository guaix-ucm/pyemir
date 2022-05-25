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
        skyref = wcsh.wcs_pix2world(pixref, 1)

    for idx, frame in enumerate(frames[1:]):
        with frame.open() as hdulist:
            wcsh = wcs.WCS(hdulist[0].header)
            pixval = wcsh.wcs_world2pix(skyref, 1)
            result[idx + 1] = -(pixval[0] - pixref[0])

    return result


def offsets_from_wcs_imgs(imgs, pixref):

    result = numpy.zeros((len(imgs), pixref.shape[1]))

    ref = imgs[0]
    wcsh = wcs.WCS(ref[0].header)
    skyref = wcsh.wcs_pix2world(pixref, 1)

    for idx, img in enumerate(imgs[1:]):
        wcsh = wcs.WCS(img[0].header)
        pixval = wcsh.wcs_world2pix(skyref, 1)
        result[idx + 1] = -(pixval[0] - pixref[0])

    return result


def reference_pix_from_wcs(frames, pixref, origin=1):
    """Compute reference pixels between frames using WCS information.

    The sky world coordinates are computed on *pixref* using
    the WCS of the first frame in the sequence. Then, the
    pixel coordinates of the reference sky world-coordinates
    are computed for the rest of the frames.

    The results is a list with the position of the reference pixel
    in each image

    """

    result = []

    with frames[0].open() as hdulist:
        wcsh = wcs.WCS(hdulist[0].header)
        skyref = wcsh.wcs_pix2world([pixref], origin)
        result.append(pixref)

    for idx, frame in enumerate(frames[1:]):
        with frame.open() as hdulist:
            wcsh = wcs.WCS(hdulist[0].header)
            pixval = wcsh.wcs_world2pix(skyref, origin)
            result.append(tuple(pixval[0]))

    return result


def reference_pix_from_wcs_imgs(imgs, pixref, origin=1):
    """Compute reference pixels between frames using WCS information.

    The sky world coordinates are computed on *pixref* using
    the WCS of the first frame in the sequence. Then, the
    pixel coordinates of the reference sky world-coordinates
    are computed for the rest of the frames.

    The results is a list with the position of the reference pixel
    in each image

    """

    result = []

    refimg = imgs[0]

    wcsh = wcs.WCS(refimg[0].header)
    skyref = wcsh.wcs_pix2world([pixref], origin)
    result.append(pixref)

    for idx, img in enumerate(imgs[1:]):
        wcsh = wcs.WCS(img[0].header)
        pixval = wcsh.wcs_world2pix(skyref, origin)
        result.append(tuple(pixval[0]))

    return result
