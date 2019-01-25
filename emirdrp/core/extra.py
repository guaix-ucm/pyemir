#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""
Extra functions
"""

import logging
import contextlib

import numpy
import astropy.io.fits as fits


from numina.array import combine


_logger = logging.getLogger(__name__)


@contextlib.contextmanager
def manage_frame(list_of_frame):
    refs = []
    for frame in list_of_frame:
        ref = frame.open()
        refs.append(ref)
    try:
        yield refs
    finally:
        # release
        for obj in refs:
            obj.close()


def combine_frames(list_of_hdulist, datamodel, method=combine.mean, errors=True, prolog=None):

    # Using the following keywords:
    # NUM-NCOM
    # UUID, TSUTC1, TSUTC2
    # VARIANCE, MAP

    import datetime
    import uuid

    cdata = list_of_hdulist

    base_header = cdata[0][0].header.copy()
    cnum = len(cdata)
    _logger.info("stacking %d images using '%s'", cnum, method.__name__)
    data = method([d[0].data for d in cdata], dtype='float32')
    hdu = fits.PrimaryHDU(data[0], header=base_header)
    _logger.debug('update result header')
    if prolog:
        _logger.debug('write prolog')
        hdu.header['history'] = prolog
    hdu.header['history'] = "Combined %d images using '%s'" % (cnum, method.__name__)
    hdu.header['history'] = 'Combination time {}'.format(datetime.datetime.utcnow().isoformat())
    for img in cdata:
        hdu.header['history'] = "Image {}".format(datamodel.get_imgid(img))
    prevnum = base_header.get('NUM-NCOM', 1)
    hdu.header['NUM-NCOM'] = prevnum * cnum
    hdu.header['UUID'] = str(uuid.uuid1())
    # Headers of last image
    hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
    if errors:
        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        result = fits.HDUList([hdu, varhdu, num])
    else:
        result = fits.HDUList([hdu])

    return result


def combine_frames2(list_of_hdulist, datamodel, method=combine.mean, errors=True, prolog=None):

    # Using the following keywords:
    # NUM-NCOM
    # UUID, TSUTC1, TSUTC2
    # VARIANCE, MAP

    import datetime
    import uuid

    cdata = list_of_hdulist

    base_header = cdata[0][0].header.copy()
    cnum = len(cdata)
    _logger.info("stacking %d images using '%s'", cnum, method.__name__)
    # Scaling, get median
    medians = numpy.array([numpy.median(d[0].data) for d in cdata])
    scales = medians / medians.max()

    data = method([d[0].data for d in cdata], dtype='float32', scales=scales)
    hdu = fits.PrimaryHDU(data[0], header=base_header)
    _logger.debug('update result header')
    if prolog:
        _logger.debug('write prolog')
        hdu.header['history'] = prolog
    hdu.header['history'] = "Combined %d images using '%s'" % (cnum, method.__name__)
    hdu.header['history'] = 'Combination time {}'.format(datetime.datetime.utcnow().isoformat())
    for img in cdata:
        hdu.header['history'] = "Image {}".format(datamodel.get_imgid(img))
    prevnum = base_header.get('NUM-NCOM', 1)
    hdu.header['NUM-NCOM'] = prevnum * cnum
    hdu.header['UUID'] = str(uuid.uuid1())
    # Headers of last image
    hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
    if errors:
        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        num = fits.ImageHDU(data[2], name='MAP')
        result = fits.HDUList([hdu, varhdu, num])
    else:
        result = fits.HDUList([hdu])

    return result


def generate_bpm_hdu(hdu):
    return generate_bpm_data(hdu.data, dtype='uint8', header=hdu.header)


def generate_empty_bpm_hdu(hdu, dtype='uint8'):
    data = numpy.zeros_like(hdu.data, dtype)
    return generate_bpm_data(data)


def generate_bpm_data(data, dtype='uint8', header=None):
    # Mask is assumed to be UINT8
    # casting here
    data = data.astype(dtype)
    hdu_bpm = fits.ImageHDU(
        data,
        header=header
    )

    hdu_bpm.header['EXTNAME'] = 'BPM'
    return hdu_bpm
