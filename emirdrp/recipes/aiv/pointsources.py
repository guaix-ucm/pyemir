#
# Copyright 2013-2016 Universidad Complutense de Madrid
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

"""Test recipe to find point sources"""

from __future__ import division

import logging

import numpy
from astropy.wcs import WCS
#import astropy.io.fits as fits
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import photutils
from numina.core import RecipeError
from numina.core import Requirement, Product, Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.core.products import ArrayType

from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
from emirdrp.products import CoordinateList2DType
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement
from .flows import basic_processing_with_combination
from .flows import init_filters_pbdfs
from .common import get_dtur_from_header


_logger = logging.getLogger('numina.recipes.emir')


class TestPointSourceRecipe(EmirRecipe):

    # Recipe Requirements
    #
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    pinhole_nominal_positions = Requirement(CoordinateList2DType,
                                            'Nominal positions of the pinholes'
                                            )
    shift_coordinates = Parameter(True, 'Use header information to'
                                  ' shift the pinhole positions from (0,0) '
                                  'to X_DTU, Y_DTU')
    box_half_size = Parameter(4, 'Half of the computation box size in pixels')
    recenter = Parameter(True, 'Recenter the pinhole coordinates')
    max_recenter_radius = Parameter(2.0, 'Maximum distance for recentering')

    # Recipe Products
    frame = Product(DataFrameType)
    positions = Product(ArrayType)
    positions_alt = Product(ArrayType)
    DTU = Product(ArrayType)
    filter = Product(str)
    readmode = Product(str)
    IPA = Product(float)
    DETPA = Product(float)
    DTUPA = Product(float)
    param_recenter = Product(bool)
    param_max_recenter_radius = Product(float)
    param_box_half_size = Product(float)

    def run(self, rinput):
        _logger.info('starting processing for slit detection')

        flow = init_filters_pbdfs(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        _logger.debug('finding point sources')

        try:
            filtername = hdr['FILTER']
            readmode = hdr['READMODE']
            ipa = hdr['IPA']
            detpa = hdr['DETPA']
            dtupa = hdr['DTUPA']
            dtub, dtur = get_dtur_from_header(hdr)
        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

        data = hdulist[0].data
        wcs = WCS(header=hdr)

        snr_detect = 5.0
        fwhm = 4.0
        npixels = 15
        box_shape = [64, 64]
        _logger.info('point source detection2')
        mask = numpy.zeros_like(data, dtype='bool')
        _logger.info('compute background map, %s', box_shape)
        bkg = photutils.background.Background(data, box_shape=box_shape,
                                              mask=mask)

        _logger.info('reference fwhm is %5.1f pixels', fwhm)
        _logger.info('detect threshold, %3.1f over background', snr_detect)
        threshold = photutils.detect_threshold(data, background=bkg.background,
                                               snr=snr_detect, mask=mask)
        _logger.info('convolve with gaussian kernel, FWHM %3.1f pixels', fwhm)
        sigma = fwhm * gaussian_fwhm_to_sigma  # FWHM = 2.

        kernel = Gaussian2DKernel(sigma)
        kernel.normalize()

        _logger.info('create segmentation image, npixels >= %d', npixels)
        segm = photutils.detect_sources(data, threshold, npixels=npixels, filter_kernel=kernel)

        # fits.writeto('segm.fits', segm.data, clobber=True)

        _logger.info('compute properties')
        props = photutils.source_properties(data - bkg.background,
                                            segm, mask=mask,
                                            background=bkg.background,
                                            wcs=wcs)

        tbl = photutils.properties_table(props)
        tbl.remove_columns(['source_sum_err',])
        # 'id', 'xcentroid', 'ycentroid', 'ra_icrs_centroid', 'dec_icrs_centroid',
        # 'source_sum', 'source_sum_err', 'background_sum', 'background_mean', 'background_at_centroid',
        # 'xmin', 'xmax', 'ymin', 'ymax', 'min_value', 'max_value', 'minval_xpos', 'minval_ypos', 'maxval_xpos', 'maxval_ypos', 'area',
        # 'equivalent_radius', 'perimeter', 'semimajor_axis_sigma', 'semiminor_axis_sigma', 'eccentricity', 'orientation', 'ellipticity',
        # 'elongation', 'covar_sigx2', 'covar_sigxy', 'covar_sigy2', 'cxx', 'cxy', 'cyy'
        # print(tbl.columns)
        positions_alt = tbl.as_array()
        positions_alt = numpy.array(positions_alt.tolist())

        result = self.create_result(frame=hdulist,
                                positions=positions_alt,
                                positions_alt=positions_alt,
                                filter=filtername,
                                DTU=dtub,
                                readmode=readmode,
                                IPA=ipa,
                                DETPA=detpa,
                                DTUPA=dtupa,
                                param_recenter=rinput.recenter,
                                param_max_recenter_radius=rinput.max_recenter_radius,
                                param_box_half_size=rinput.box_half_size
                                )
        return result
