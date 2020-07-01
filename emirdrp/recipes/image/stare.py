#
# Copyright 2011-2020 Universidad Complutense de Madrid
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


from numina.array import combine
from numina.core import Result
from numina.core.query import Ignore
from numina.core.recipes import timeit
from numina.processing.combine import basic_processing_with_combination
from numina.util.context import manage_fits

from emirdrp.core.recipe import EmirRecipe
import emirdrp.core.extra as extra
import emirdrp.requirements as reqs
import emirdrp.products as prods
import emirdrp.processing.combine as comb
import emirdrp.decorators


_logger = logging.getLogger(__name__)


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

        with manage_fits(frames) as list_of:
            c_img = comb.combine_images(list_of, method=combine.mean, errors=False)

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
            hdu_bpm = extra.generate_bpm_hdu(hdul_bpm[0])
        else:
            self.logger.debug('using empty BPM')
            hdu_bpm = extra.generate_empty_bpm_hdu(processed_img[0])

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
            hdu_bpm = extra.generate_bpm_hdu(hdul_bpm[0])
        else:
            hdu_bpm = extra.generate_empty_bpm_hdu(hdulist[0])

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
