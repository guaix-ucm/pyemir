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
import sep

import astropy.io.fits as fits
import numpy
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from numina.array.fwhm import compute_fwhm_2d_simple
from numina.core import Product, Parameter
from numina.core import RecipeError
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from scipy.spatial import KDTree

import emirdrp.datamodel as datamodel
from emirdrp.core import EmirRecipe
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.products import DataFrameType
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement
from .procedures import image_box2d

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
    ROTANG = Product(float)
    DETPA = Product(float)
    DTUPA = Product(float)
    param_recenter = Product(bool)
    param_max_recenter_radius = Product(float)
    param_box_half_size = Product(float)

    def run(self, rinput):
        _logger.info('starting processing for object detection')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        _logger.debug('finding point sources')

        try:
            filtername = hdr['FILTER']
            readmode = hdr['READMODE']
            rotang = hdr['ROTANG']
            detpa = hdr['DETPA']
            dtupa = hdr['DTUPA']
            dtub, dtur = datamodel.get_dtur_from_header(hdr)
        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

        data = hdulist[0].data

        # Copy needed in numpy 1.7
        # This seems already bitswapped??
        # FIXME: check this works offline/online
        # ndata = data.byteswap().newbyteorder()
        # data = data.byteswap(inplace=True).newbyteorder()

        snr_detect = 5.0
        fwhm = 4.0
        npixels = 15
        box_shape = [64, 64]
        _logger.info('point source detection2')
        _logger.info('using internal mask to remove corners')
        # Corners
        mask = numpy.zeros_like(data, dtype='int32')
        mask[2000:, 0:80] = 1
        mask[2028:, 2000:] = 1
        mask[:50, 1950:] = 1
        mask[:100, :50] = 1
        # Remove corner regions

        _logger.info('compute background map, %s', box_shape)
        bkg = sep.Background(data)

        _logger.info('reference fwhm is %5.1f pixels', fwhm)
        _logger.info('detect threshold, %3.1f over background', snr_detect)
        _logger.info('convolve with gaussian kernel, FWHM %3.1f pixels', fwhm)
        sigma = fwhm * gaussian_fwhm_to_sigma
        #
        kernel = Gaussian2DKernel(sigma)
        kernel.normalize()

        thresh = snr_detect * bkg.globalrms
        data_s = data - bkg.back()
        objects, segmap = sep.extract(data - bkg.back(), thresh, minarea=npixels,
                                      filter_kernel=kernel.array, segmentation_map=True,
                                      mask=mask)
        fits.writeto('segmap.fits', segmap)
        _logger.info('detected %d objects', len(objects))

        # Hardcoded values
        rs2 = 15.0
        fit_rad = 10.0
        flux_min = 1000.0
        flux_max = 30000.0
        _logger.debug('Flux limit is %6.1f %6.1f', flux_min, flux_max)
        # FIXME: this should be a view, not a copy
        xall = objects['x']
        yall = objects['y']
        mm = numpy.array([xall, yall]).T
        _logger.info('computing FWHM')
        # Find objects with pairs inside fit_rad
        kdtree = KDTree(mm)
        nearobjs = (kdtree.query_ball_tree(kdtree, r=fit_rad))
        positions = []
        for idx, obj in enumerate(objects):
            x0 = obj['x']
            y0 = obj['y']
            sl = image_box2d(x0, y0, data.shape, (fit_rad, fit_rad))
            # sl_sky = image_box2d(x0, y0, data.shape, (rs2, rs2))
            part_s = data_s[sl]
            # Logical coordinates
            xx0 = x0 - sl[1].start
            yy0 = y0 - sl[0].start

            _, fwhm_x, fwhm_y = compute_fwhm_2d_simple(part_s, xx0, yy0)

            if min(fwhm_x, fwhm_x) < 3:
                continue
            if flux_min > obj['peak'] or flux_max < obj['peak']:
                continue
            # nobjs is the number of object inside fit_rad
            nobjs = len(nearobjs[idx])

            flag = 0 if nobjs == 1 else 1

            positions.append([idx, x0, y0, obj['peak'], fwhm_x, fwhm_y, flag])

        _logger.info('saving photometry')
        positions = numpy.array(positions)
        positions_alt = positions
        _logger.info('end processing for object detection')

        result = self.create_result(frame=hdulist,
                                positions=positions_alt,
                                positions_alt=positions_alt,
                                filter=filtername,
                                DTU=dtub,
                                readmode=readmode,
                                ROTANG=rotang,
                                DETPA=detpa,
                                DTUPA=dtupa,
                                param_recenter=rinput.recenter,
                                param_max_recenter_radius=rinput.max_recenter_radius,
                                param_box_half_size=rinput.box_half_size
                                )
        return result

