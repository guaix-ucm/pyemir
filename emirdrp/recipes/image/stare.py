#
# Copyright 2011-2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


"""
Stare Image mode of EMIR
"""

import logging
import contextlib

import numpy
import astropy.io.fits as fits


from numina.array import combine
from numina.core import Result
from numina.core.query import Ignore
from numina.core.recipes import timeit

from emirdrp.core.recipe import EmirRecipe
import emirdrp.requirements as reqs
import emirdrp.products as prods
from emirdrp.processing.combine import basic_processing_with_combination
import emirdrp.decorators


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


class StareImageRecipe2(EmirRecipe):
    """Process images in Stare Image Mode"""

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()

    result_image = Result(prods.ProcessedImage)

    def run(self, rinput):
        self.logger.info('starting stare image reduction (offline)')

        frames = rinput.obresult.frames
        datamodel = self.datamodel

        with manage_frame(frames) as list_of:
            c_img = combine_frames(list_of, datamodel, method=combine.mean, errors=False)

        flow = self.init_filters(rinput)

        # Correct Bias if needed
        # Correct Dark if needed
        # Correct FF

        processed_img = flow(c_img)

        hdr = processed_img[0].header
        self.set_base_headers(hdr)

        self.logger.debug('append BPM')
        if rinput.master_bpm is not None:
            self.logger.debug('using BPM from inputs')
            hdul_bpm = rinput.master_bpm.open()
            hdu_bpm = generate_bpm_hdu(hdul_bpm[0])
        else:
            self.logger.debug('using empty BPM')
            hdu_bpm = generate_empty_bpm_hdu(processed_img[0])

        # Append the BPM to the result
        processed_img.append(hdu_bpm)
        self.logger.info('end stare image (off) reduction')
        result = self.create_result(result_image=processed_img)

        return result


    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(StareImageRecipe2, self).set_base_headers(hdr)
        # Set EXP to 0
        hdr['EXP'] = 0
        return hdr


class StareImageBaseRecipe(EmirRecipe):
    """Process images in Stare Image Mode"""

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement(optional=True)

    frame = Result(prods.ProcessedImage)

    def __init__(self, *args, **kwargs):
        super(StareImageBaseRecipe, self).__init__(*args, **kwargs)
        if False:
            self.query_options['master_sky'] = Ignore()

    @emirdrp.decorators.loginfo
    @timeit
    def run(self, rinput):
        self.logger.info('starting stare image reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(
            rinput,
            flow,
            method=combine.median
        )
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        if rinput.master_bpm:
            hdul_bpm = rinput.master_bpm.open()
            hdu_bpm = generate_bpm_hdu(hdul_bpm[0])
        else:
            hdu_bpm = generate_empty_bpm_hdu(hdulist[0])

        # Append the BPM to the result
        hdulist.append(hdu_bpm)
        self.logger.info('end stare image reduction')
        result = self.create_result(frame=hdulist)

        return result

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(StareImageBaseRecipe, self).set_base_headers(hdr)
        # Update EXP to 0
        hdr['EXP'] = 0
        return hdr


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
