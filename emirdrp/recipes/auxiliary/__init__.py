#
# Copyright 2011-2016 Universidad Complutense de Madrid
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

"""Auxiliary Recipes for EMIR"""

import logging

import numpy
from numina.array.combine import median
from numina.core import DataFrame
from numina.core import Product
from numina.core import RecipeError
from numina.core.requirements import ObservationResultRequirement

import emirdrp.instrument.channels as allchannels
from emirdrp.core import EMIR_BIAS_MODES
from emirdrp.core import EmirRecipe
from emirdrp.processing.info import gather_info_frames
from emirdrp.products import ChannelLevelStatistics
from emirdrp.products import ChannelLevelStatisticsType
from emirdrp.products import MasterBias, MasterDark
from emirdrp.products import MasterIntensityFlat
from emirdrp.products import SlitTransmissionCalibration
from emirdrp.products import WavelengthCalibration, MasterSpectralFlat
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSpectralFlatFieldRequirement
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.processing.combine import  basic_processing_with_segmentation


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
    
    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()

    biasframe = Product(MasterBias)

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

    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()

    darkframe = Product(MasterDark)

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

    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    flatframe = Product(MasterIntensityFlat)    

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


class SimpleSkyRecipe(EmirRecipe):
    """Recipe to process data taken in intensity flat-field mode.

    """

    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()

    skyframe = Product(MasterIntensityFlat)


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

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()

    skyframe = Product(MasterIntensityFlat)


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

    master_bpm = MasterBadPixelMaskRequirement()
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    
    flatframe = Product(MasterSpectralFlat)

    def run(self, rinput):
        return self.create_result(flatframe=MasterSpectralFlat())


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

    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    slit = Product(SlitTransmissionCalibration)

    def run(self, rinput):
        return self.create_result(slit=SlitTransmissionCalibration())



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

    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_spectral_ff = MasterSpectralFlatFieldRequirement()

    cal = Product(WavelengthCalibration)

    def run(self, rinput):
        return self.create_result(cal=WavelengthCalibration())
