#
# Copyright 2013-2014 Universidad Complutense de Madrid
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

'''Recipe to detect slits in the AIV mask'''

from __future__ import division
#
import logging
# import math
#
# import numpy
# import scipy.interpolate as itpl
# import scipy.optimize as opz
# from astropy.modeling import models, fitting
from astropy.io import fits
# import photutils
#
# from numina.array.recenter import img_box, centering_centroid
#
from numina import __version__
from numina.core import BaseRecipe, RecipeRequirements, RecipeError
from numina.core import Requirement, Product, DataProductRequirement, Parameter
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement
from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import FlatFieldCorrector, SkyCorrector
from numina.flow import SerialFlow
from numina.flow.node import IdNode
from numina.array import combine
#
from emir.core import RecipeResult
from emir.core import EMIR_BIAS_MODES
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import DataFrameType, MasterIntensityFlat
# from emir.dataproducts import CoordinateList2DType
# from emir.dataproducts import ArrayType
from emir.core import gather_info
# from photutils import aperture_circular
#
# from .procedures import compute_fwhm_spline2d_fit
# from .procedures import compute_fwhm_enclosed_direct
# from .procedures import compute_fwhm_enclosed_grow
# from .procedures import compute_fwhm_simple
# from .procedures import moments
# from .procedures import AnnulusBackgroundEstimator
# from .procedures import img_box2d

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"

GAUSS_FWHM_FACTOR = 2.354800


class TestSlitDetectionRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')
    master_sky = DataProductRequirement(MasterIntensityFlat, 'Master Sky calibration')
#    pinhole_nominal_positions = Requirement(CoordinateList2DType, 'Nominal positions of the pinholes')
#    shift_coordinates = Parameter(True, 'Use header information to shift the pinhole positions from (0,0) to X_DTU, Y_DTU')
#    box_half_size = Parameter(4, 'Half of the computation box size in pixels')
#    recenter = Parameter(True, 'Recenter the pinhole coordinates')
#    max_recenter_radius = Parameter(2.0, 'Maximum distance for recentering')


class TestSlitDetectionRecipeResult(RecipeResult):
    frame = Product(DataFrameType)
#    positions = Product(ArrayType)
#    positions_alt = Product(ArrayType)
#    DTU = Product(ArrayType)
#    filter = Product(str)
#    readmode = Product(str)
#    IPA = Product(float)


@define_requirements(TestSlitDetectionRecipeRequirements)
@define_result(TestSlitDetectionRecipeResult)
class TestSlitDetectionRecipe(BaseRecipe):

    def __init__(self):
        super(TestSlitDetectionRecipe, self).__init__(
            author=_s_author,
            version="0.1.0"
            )

    def run(self, rinput):
        _logger.info('starting pinhole processing')

        meta = gather_info(rinput)
        iinfo = meta['obresult']

        if iinfo:
            mode = iinfo[0]['readmode']
            if mode.lower() in EMIR_BIAS_MODES:
                use_bias = True
                _logger.info('readmode is %s, bias required', mode)

            else:
                use_bias = False
                _logger.info('readmode is %s, no bias required', mode)

        dark_info = meta['master_dark']
        flat_info = meta['master_flat']
        sky_info = meta['master_sky']

        print('images info:', iinfo)
        if use_bias:
            bias_info = meta['master_bias']
            print('bias info:', bias_info)
        print('dark info:', dark_info)
        print('flat info:', flat_info)
        print('sky info:', sky_info)

        # Loading calibrations
        if use_bias:
            with rinput.master_bias.open() as hdul:
                _logger.info('loading bias')
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbias)
        else:
            _logger.info('ignoring bias')
            bias_corrector = IdNode()

        with rinput.master_dark.open() as mdark_hdul:
            _logger.info('loading dark')
            mdark = mdark_hdul[0].data
            dark_corrector = DarkCorrector(mdark)

        with rinput.master_flat.open() as mflat_hdul:
            _logger.info('loading intensity flat')
            mflat = mflat_hdul[0].data
            flat_corrector = FlatFieldCorrector(mflat)

        with rinput.master_sky.open() as msky_hdul:
            _logger.info('loading sky')
            msky = msky_hdul[0].data
            sky_corrector = SkyCorrector(msky)

        flow = SerialFlow([bias_corrector, dark_corrector,
                           flat_corrector, sky_corrector])

        odata = []
        cdata = []
        try:
            _logger.info('processing input frames')
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                fname = hdulist.filename()
                if fname:
                    _logger.info('input is %s', fname)
                else:
                    _logger.info('input is %s', hdulist)

                final = flow(hdulist)
                _logger.debug('output is input: %s', final is hdulist)

                cdata.append(final)

                # Files to be closed at the end
                odata.append(hdulist)
                if final is not hdulist:
                    odata.append(final)

            _logger.info("stacking %d images using 'mean'", len(cdata))
            data = combine.mean([d[0].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0][0].header.copy())

        finally:
            _logger.debug('closing images')
            for hdulist in odata:
                hdulist.close()

        _logger.debug('update result header')
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')

        _logger.debug('finding pinholes')

        try:
            filtername = hdr['FILTER']
            readmode = hdr['READMODE']
            ipa = hdr['IPA']
            xdtu = hdr['XDTU']
            ydtu = hdr['YDTU']
            zdtu = hdr['ZDTU']
        except KeyError as error:
            _logger.error(error)
            raise RecipeError(error)

        result = self.create_result(frame=hdulist)

        return result
