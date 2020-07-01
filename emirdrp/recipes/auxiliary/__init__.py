#
# Copyright 2011-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Auxiliary Recipes for EMIR"""

import logging

from numina.array.combine import median
from numina.core import DataFrame
from numina.core import Result
from numina.exceptions import RecipeError
from numina.util.context import manage_fits

from emirdrp.core import EMIR_BIAS_MODES
from emirdrp.core.recipe import EmirRecipe
from emirdrp.processing.info import gather_info_frames
import emirdrp.products as prods
import emirdrp.requirements as reqs

from numina.processing.combine import basic_processing_with_combination
from emirdrp.processing.combine import basic_processing_with_segmentation
from emirdrp.processing.combine import combine_images, scale_with_median


_logger = logging.getLogger('numina.recipes.emir')


class BiasRecipe(EmirRecipe):
    """
    Recipe to process data taken in Bias image Mode.

    Bias images only appear in Simple Readout mode.

    **Observing modes:**

     * Bias Image (3.1)

    **Inputs:**

    **Outputs:**

     * A combined bias frame, with variance extension.
     * Statistics of the final image per channel (mean, median, variance)

    **Procedure:**

    The list of images can be readly processed by
    combining them with a median algorithm.

    """
    
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    obresult = reqs.ObservationResultRequirement()

    biasframe = Result(prods.MasterBias)

    def run(self, rinput):
        _logger.info('starting bias reduction')

        iinfo = gather_info_frames(rinput.obresult.frames)

        if iinfo:
            mode = iinfo[0]['readmode']
            if mode.lower() not in EMIR_BIAS_MODES:
                msg = 'readmode %s, is not a bias mode' % mode
                _logger.error(msg)
                raise RecipeError(msg)

        flow = lambda x: x
        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=False)

        pdata = hdulist[0].data

        # update hdu header with
        # reduction keywords
        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        hdr['CCDMEAN'] = pdata.mean()

        _logger.info('bias reduction ended')

        result = self.create_result(biasframe=DataFrame(hdulist))
        return result


class DarkRecipe(EmirRecipe):
    """Recipe to process data taken in Dark current image Mode.

    Recipe to process dark images. The dark images will be combined
    using the median.
    They do have to be of the same exposure time t.

    **Observing mode:**

     * Dark current Image (3.2)

    **Inputs:**

    **Outputs:**

     * A combined dark frame, with variance extension.
    """

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    obresult = reqs.ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()

    darkframe = Result(prods.MasterDark)

    def run(self, rinput):

        _logger.info('starting dark reduction')

        flow = self.init_filters(rinput)

        iinfo = gather_info_frames(rinput.obresult.frames)
        ref_exptime = 0.0
        for el in iinfo[1:]:
            if abs(el['texp'] - ref_exptime) > 1e-4:
                _logger.error('image with wrong exposure time')
                raise RecipeError('image with wrong exposure time')

        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=True)

        pdata = hdulist[0].data

        # update hdu header with
        # reduction keywords

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        hdr['CCDMEAN'] = pdata.mean()

        _logger.info('dark reduction ended')
        result = self.create_result(darkframe=hdulist)
        return result


class IntensityFlatRecipe(EmirRecipe):
    """Recipe to process data taken in intensity flat-field mode.

    Recipe to process intensity flat-fields. The flat-on and
    flat-off images are combined (method?) separately and the subtracted
    to obtain a thermal subtracted flat-field.

    **Observing modes:**

     * Intensity Flat-Field

    **Inputs:**

      * A master dark frame
      * Non linearity
      * A model of the detector.

    **Outputs:**

     * TBD

    **Procedure:**

     * A combined thermal subtracted flat field, normalized to median 1,
       with with variance extension and quality flag.

    """

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    obresult = reqs.ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    flatframe = Result(prods.MasterIntensityFlat)

    def run(self, rinput):
        _logger.info('starting flat reduction')

        errors = True

        flow = self.init_filters(rinput)
        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=errors)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        mm = hdulist[0].data.mean()
        hdr['CCDMEAN'] = mm

        hdulist[0].data /= mm
        if errors:
            hdulist['variance'].data /= (mm*mm)

        result = self.create_result(flatframe=hdulist)

        return result


class IntensityFlatRecipe2(EmirRecipe):
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    master_flatframe = Result(prods.MasterIntensityFlat)

    def run(self, rinput):
        from numina.array import combine

        _logger.info('starting flat reduction')

        frames = rinput.obresult.frames

        with manage_fits(frames) as list_of:
            scaled_mean = scale_with_median(combine.mean)
            c_img = combine_images(list_of, method=scaled_mean, errors=False)

        self.save_intermediate_img(c_img, 'p0.fits')

        flow = self.init_filters(rinput)

        processed_img = flow(c_img)

        self.save_intermediate_img(processed_img, 'p1.fits')

        hdr = processed_img[0].header
        self.set_base_headers(hdr)

        import scipy.ndimage.filters

        _logger.info('median filter')
        data_smooth = scipy.ndimage.filters.median_filter(processed_img[0].data, size=11)

        self.save_intermediate_array(data_smooth, 'smooth.fits')

        mm = processed_img[0].data.mean()
        hdr['CCDMEAN'] = mm

        processed_img[0].data /= data_smooth

        self.save_intermediate_img(processed_img, 'p2.fits')

        result = self.create_result(master_flatframe=processed_img)

        return result


class SimpleSkyRecipe(EmirRecipe):
    """Recipe to process data taken in intensity flat-field mode.

    """

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    obresult = reqs.ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()

    skyframe = Result(prods.MasterSky)


    def run(self, rinput):
        _logger.info('starting sky reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow,
                                                    method=median,
                                                    errors=True)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        result = self.create_result(skyframe=hdulist)

        return result


class DitherSkyRecipe(EmirRecipe):
    """Recipe to process data taken in dither sky mode.

    """

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()

    skyframe = Result(prods.MasterSky)


    def run(self, rinput):
        _logger.debug('instrument %s, mode %s', rinput.obresult.instrument,
                      rinput.obresult.mode
                      )
        _logger.info('starting sky reduction with dither')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_segmentation(rinput, flow,
                                                    method=median,
                                                    errors=True)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        _logger.info('end sky reduction with dither')

        result = self.create_result(skyframe=hdulist)

        return result


class SpectralFlatRecipe(EmirRecipe):
    """Recipe to process data taken in intensity flat-field mode."""

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    obresult = reqs.ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    
    flatframe = Result(prods.MasterSpectralFlat)

    def run(self, rinput):
        return self.create_result(flatframe=prods.MasterSpectralFlat())


class SlitTransmissionRecipe(EmirRecipe):
    """Recipe to calibrate the slit transmission.

    **Observing modes:**

        * Slit transmission calibration (4.4)

    **Inputs:**

        * A list of uniformly illuminated images of MSM

    **Outputs:**

     * A list of slit transmission functions

    **Procedure:**

     * TBD

    """

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()

    slit = Result(prods.SlitTransmissionCalibration)

    def run(self, rinput):
        return self.create_result(slit=prods.SlitTransmissionCalibration())



class WavelengthCalibrationRecipe(EmirRecipe):
    """Recipe to calibrate the spectral response.

    **Observing modes:**

        * Wavelength calibration (4.5)

    **Inputs:**

     * List of line positions
     * Calibrations up to spectral flatfielding

    **Outputs:**

     * Wavelength calibration structure

    **Procedure:**

     * TBD
    """

    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_spectral_ff = reqs.MasterSpectralFlatFieldRequirement()

    cal = Result(prods.WavelengthCalibration)

    def run(self, rinput):
        return self.create_result(cal=prods.WavelengthCalibration())
