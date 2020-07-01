#
# Copyright 2013-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Combination routines"""

from __future__ import division

import datetime
import functools
import logging
import uuid

import numpy
from astropy.io import fits
import numina.datamodel as datam
from numina.frame.utils import copy_img
from numina.array import combine
from numina.array import combine_shape

from emirdrp.processing.wcs import offsets_from_wcs


_logger = logging.getLogger(__name__)


def basic_processing(rinput, flow):
    """

    Parameters
    ----------
    rinput : RecipeInput
    flow

    Returns
    -------

    """
    return basic_processing_(rinput.obresult.images, flow)


def process_abba(images, errors=False, prolog=None):
    """
    Process images in ABBA sequence

    Parameters
    ----------
    images
    errors
    prolog

    Returns
    -------

    """
    cnum = len(images)
    result = copy_img(images[0])
    base_header = result[0].header.copy()

    now = datetime.datetime.utcnow()
    _logger.info("processing ABBA")
    ###
    ###
    dataA0 = images[0][0].data.astype('float32')
    dataB0 = images[1][0].data.astype('float32')
    dataB1 = images[2][0].data.astype('float32')
    dataA1 = images[3][0].data.astype('float32')
    dataAB0 = dataA0 - dataB0
    dataAB1 = dataA1 - dataB1
    data = dataAB0 + dataAB1
    ###
    hdu = result[0]
    hdu.data = data
    _logger.debug('update result header')
    if prolog:
        hdu.header['history'] = prolog
    hdu.header['history'] = "Processed ABBA"
    hdu.header['history'] = 'Combination time {}'.format(now.isoformat())
    for img in images:
        hdu.header['history'] = "Image {}".format(datam.get_imgid(img))
    prevnum = base_header.get('NUM-NCOM', 1)
    hdu.header['NUM-NCOM'] = prevnum * cnum // 2
    hdu.header['UUID'] = str(uuid.uuid1())
    # Headers of last image
    hdu.header['TSUTC2'] = images[-1][0].header['TSUTC2']
    # Not sure why this is needed
    hdu.header['EXTEND'] = True
    for img, key in zip(images, ['A', 'B', 'B', 'A']):
        imgid = datam.get_imgid(img)
        hdu.header['history'] = "Image '{}' is '{}'".format(imgid, key)

    # Update primary extension
    result[0] = hdu

    return result


def process_ab(images, errors=False, prolog=None):
    """
    Process images in AB sequence

    Parameters
    ----------
    images
    errors
    prolog

    Returns
    -------

    """
    cnum = len(images)
    result = copy_img(images[0])
    base_header = result[0].header.copy()
    now = datetime.datetime.utcnow()
    _logger.info("processing AB")
    ###
    ###
    dataA0 = images[0][0].data.astype('float32')
    dataB0 = images[1][0].data.astype('float32')
    data = dataA0 - dataB0

    ###
    hdu = fits.PrimaryHDU(data, header=base_header)
    _logger.debug('update result header')
    if prolog:
        hdu.header['history'] = prolog
    hdu.header['history'] = "Processed AB"
    hdu.header['history'] = 'Combination time {}'.format(now.isoformat())
    for img in images:
        hdu.header['history'] = "Image {}".format(datam.get_imgid(img))
    prevnum = base_header.get('NUM-NCOM', 1)
    hdu.header['NUM-NCOM'] = prevnum * cnum // 2
    hdu.header['UUID'] = str(uuid.uuid1())
    # Headers of last image
    hdu.header['TSUTC2'] = images[-1][0].header['TSUTC2']

    for img, key in zip(images, ['A', 'B']):
        imgid = datam.get_imgid(img)
        hdu.header['history'] = "Image '{}' is '{}'".format(imgid, key)

    # Update primary extension
    # Not sure why this is needed
    hdu.header['EXTEND'] = True
    result[0] = hdu

    return result


def create_proc_hdulist(images, data_array):
    cnum = len(images)
    result = copy_img(images[0])
    base_header = result[0].header.copy()

    hdu = fits.PrimaryHDU(data_array, header=base_header)
    # self.set_base_headers(hdu.header)
    hdu.header['UUID'] = str(uuid.uuid1())
    # Update obsmode in header
    #hdu.header['OBSMODE'] = 'LS_ABBA'
    # Headers of last image
    hdu.header['TSUTC2'] = images[-1][0].header['TSUTC2']
    # Not sure why this is needed
    hdu.header['EXTEND'] = True
    result[0] = hdu
    return result


def basic_processing_(frames, flow):
    """

    Parameters
    ----------
    frames
    flow

    Returns
    -------

    """

    cdata = []
    _logger.info('processing input frames')
    for frame in frames:
        hdulist = frame.open()
        fname = datam.get_imgid(hdulist)
        _logger.info('input is %s', fname)
        final = flow(hdulist)
        _logger.debug('output is input: %s', final is hdulist)

        cdata.append(final)

    return cdata


def scale_with_median(method):
    """Adapt combine method to scale with median
    """
    @functools.wraps(method)
    def wrapper(data_list, dtype='float32'):
        medians = numpy.array([numpy.median(d) for d in data_list])
        scales = medians / medians.max()
        data_com = method(data_list, dtype=dtype, scales=scales)
        return data_com
    # Adapt __name__
    rename = "scaled_" + wrapper.__name__
    wrapper.__name__ = rename

    return wrapper


# FIXME: use the function in numina
def combine_images(images, method=combine.mean, method_kwargs=None,
                   errors=False, prolog=None):
    """Combine a sequence of HDUList objects.

    Using the following keywords:
     * NUM-NCOM
     * UUID, TSUTC1, TSUTC2
     * VARIANCE, MAP

    Parameters
    ----------
    images
    method
    method_kwargs : dict (optional)
    errors : bool (optional)
    prolog

    Returns
    -------

    """
    cnum = len(images)
    result = copy_img(images[0])
    base_header = result[0].header.copy()

    now = datetime.datetime.utcnow()
    _logger.info("stacking %d images using '%s'", cnum, method.__name__)
    data = method([d[0].data for d in images], dtype='float32')
    hdu = result[0]
    hdu.data = data[0]
    _logger.debug('update result header')
    if prolog:
        hdu.header['history'] = prolog
    hdu.header['history'] = "Combined %d images using '%s'" % (cnum, method.__name__)
    hdu.header['history'] = 'Combination time {}'.format(now.isoformat())
    for img in images:
        hdu.header['history'] = "Image {}".format(datam.get_imgid(img))
    prevnum = base_header.get('NUM-NCOM', 1)
    hdu.header['NUM-NCOM'] = prevnum * cnum
    hdu.header['UUID'] = str(uuid.uuid1())
    # Headers of last image
    hdu.header['TSUTC2'] = images[-1][0].header['TSUTC2']
    # Not sure why this is needed
    hdu.header['EXTEND'] = True
    result[0] = hdu
    if errors:
        varhdu = fits.ImageHDU(data[1], name='VARIANCE')
        result.append(varhdu)
        num = fits.ImageHDU(data[2].astype('int16'), name='MAP')
        result.append(num)
    return result


def resize_hdul(hdul, newshape, region, extensions=None, window=None,
                    scale=1, fill=0.0, conserve=True):
    """

    Parameters
    ----------
    hdul
    newshape
    region
    extensions
    window
    scale
    fill
    conserve

    Returns
    -------

    """
    from numina.frame import resize_hdu

    if extensions is None:
        extensions = [0]

    nhdul = [None] * len(hdul)
    for ext, hdu in enumerate(hdul):
        if ext in extensions:
            nhdul[ext] = resize_hdu(hdu, newshape,
                                    region, fill=fill,
                                    window=window,
                                    scale=scale,
                                    conserve=conserve)
        else:
            nhdul[ext] = hdu
    return fits.HDUList(nhdul)


def resize(frames, shape, offsetsp, finalshape, window=None):
    """

    Parameters
    ----------
    frames
    shape
    offsetsp
    finalshape
    window

    Returns
    -------

    """
    from numina.array import subarray_match
    _logger.info('Resizing frames and masks')
    rframes = []
    regions = []
    for frame, rel_offset in zip(frames, offsetsp):
        region, _ = subarray_match(finalshape, rel_offset, shape)
        rframe = resize_hdul(frame.open(), finalshape, region)
        rframes.append(rframe)
        regions.append(region)
    return rframes, regions


def resize_hdulists(hdulists, shape, offsetsp, finalshape, window=None):
    """

    Parameters
    ----------
    hdulists
    shape
    offsetsp
    finalshape
    window

    Returns
    -------

    """
    from numina.array import subarray_match
    _logger.info('Resizing frames and masks')
    rhdulist = []
    regions = []
    for hdulist, rel_offset in zip(hdulists, offsetsp):
        region, _ = subarray_match(finalshape, rel_offset, shape)
        rframe = resize_hdul(hdulist, finalshape, region)
        rhdulist.append(rframe)
        regions.append(region)
    return rhdulist, regions


def basic_processing_with_segmentation(rinput, flow,
                                          method=combine.mean,
                                          errors=True, bpm=None):

    odata = []
    cdata = []
    try:
        _logger.info('processing input images')
        for frame in rinput.obresult.images:
            hdulist = frame.open()
            fname = datam.get_imgid(hdulist)
            _logger.info('input is %s', fname)
            final = flow(hdulist)
            _logger.debug('output is input: %s', final is hdulist)

            cdata.append(final)

            # Files to be closed at the end
            odata.append(hdulist)
            if final is not hdulist:
                odata.append(final)

        base_header = cdata[0][0].header.copy()

        baseshape = cdata[0][0].data.shape
        subpixshape = cdata[0][0].data.shape

        _logger.info('Computing offsets from WCS information')
        refpix = numpy.divide(numpy.array([baseshape], dtype='int'), 2).astype('float')
        offsets_xy = offsets_from_wcs(rinput.obresult.frames, refpix)
        _logger.debug("offsets_xy %s", offsets_xy)
        # Offsets in numpy order, swaping
        offsets_fc = offsets_xy[:, ::-1]
        offsets_fc_t = numpy.round(offsets_fc).astype('int')

        _logger.info('Computing relative offsets')
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)
        _logger.debug("offsetsp %s", offsetsp)

        _logger.info('Shape of resized array is %s', finalshape)
        # Resizing target frames
        rhduls, regions = resize_hdulists(cdata, subpixshape, offsetsp, finalshape)

        _logger.info("stacking %d images, with offsets using '%s'", len(cdata), method.__name__)
        data1 = method([d[0].data for d in rhduls], dtype='float32')

        segmap = segmentation_combined(data1[0])
        # submasks
        if bpm is None:
            masks = [(segmap[region] > 0) for region in regions]
        else:
            masks = [((segmap[region] > 0) & bpm) for region in regions]

        _logger.info("stacking %d images, with objects mask using '%s'", len(cdata), method.__name__)
        data2 = method([d[0].data for d in cdata], masks=masks, dtype='float32')
        hdu = fits.PrimaryHDU(data2[0], header=base_header)
        points_no_data = (data2[2] == 0).sum()

        _logger.debug('update result header')
        hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
        hdu.header['history'] = "Combined %d images using '%s'" % (len(cdata), method.__name__)
        hdu.header['history'] = 'Combination time {}'.format(datetime.datetime.utcnow().isoformat())
        hdu.header['UUID'] = str(uuid.uuid1())
        _logger.info("missing points, total: %d, fraction: %3.1f", points_no_data, points_no_data / data2[2].size)

        if errors:
            varhdu = fits.ImageHDU(data2[1], name='VARIANCE')
            num = fits.ImageHDU(data2[2], name='MAP')
            result = fits.HDUList([hdu, varhdu, num])
        else:
            result = fits.HDUList([hdu])
    finally:
        _logger.debug('closing images')
        for hdulist in odata:
            hdulist.close()

    return result


def segmentation_combined(data, snr_detect=10.0, fwhm=4.0, npixels=15, mask_corners=False):
    import sep
    from astropy.convolution import Gaussian2DKernel
    from astropy.stats import gaussian_fwhm_to_sigma

    box_shape = [64, 64]
    _logger.info('point source detection')

    # Corners
    mask = numpy.zeros_like(data, dtype='int32')
    if mask_corners:
        # Remove corner regions
        _logger.info('using internal mask to remove corners')
        mask[2000:, 0:80] = 1
        mask[2028:, 2000:] = 1
        mask[:50, 1950:] = 1
        mask[:100, :50] = 1

    _logger.info('compute background map, %s', box_shape)
    bkg = sep.Background(data)

    _logger.info('reference fwhm is %5.1f pixels', fwhm)

    _logger.info('convolve with gaussian kernel, FWHM %3.1f pixels', fwhm)
    sigma = fwhm * gaussian_fwhm_to_sigma
    #
    kernel = Gaussian2DKernel(sigma)
    kernel.normalize()

    _logger.info('background level is %5.1f', bkg.globalback)
    _logger.info('background rms is %5.1f', bkg.globalrms)
    _logger.info('detect threshold, %3.1f sigma over background', snr_detect)
    thresh = snr_detect * bkg.globalrms

    data_s = data - bkg.back()
    try:
        objects, segmap = sep.extract(data_s, thresh, minarea=npixels,
                                      filter_kernel=kernel.array, segmentation_map=True,
                                      mask=mask)
        _logger.info('detected %d objects', len(objects))
    except Exception as error:
        _logger.warning("%s", error)
        segmap = numpy.zeros_like(data_s, dtype='int')
    return segmap


def basic_processing_with_update(rinput, flow):

    # FIXME: this only works with local images
    # We don't know how to store temporary GCS frames
    _logger.info('processing input images')
    for frame in rinput.obresult.images:
        with fits.open(frame.label, mode='update') as hdul:
            flow(hdul)
