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

from __future__ import division, print_function

import logging

import numpy
from numina.array.peaks.findpeaks1D import findPeaks_spectrum, refinePeaks_spectrum
from numina.array.wavecal.arccalibration import arccalibration_direct, fit_solution, \
                                        gen_triplets_master
from numina.array.wavecal.statsummary import sigmaG
from numina.core import Requirement, Product, Parameter
from numina.core.products import ArrayType
from numina.core.products import LinesCatalog
from numina.core.requirements import ObservationResultRequirement
from scipy.interpolate import interp1d

from emirdrp.core import EmirRecipe
from emirdrp.products import SlitsCatalog
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.processing.combine import basic_processing_with_combination


_logger = logging.getLogger('numina.recipes.emir')


class ArcCalibrationRecipe(EmirRecipe):

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    lines_catalog = Requirement(LinesCatalog, "Catalog of lines")
    slits_catalog= Requirement(SlitsCatalog, "Catalog of slits")
    polynomial_degree = Parameter(2, 'Polynomial degree of the arc calibration')

    polynomial_coeffs = Product(ArrayType)

    def run(self, rinput):
        _logger.info('starting arc calibration')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)


        # ToDo
        # - replace TBM (To Be Modified) by proper code
        # - eliminate TBR code (To Be Removed)
        # - Set LDEBUG and LPLOT to False
        # - read slitletdef.txt to read slitlet definition for current image
        # - loop in nslits to excecute the code in calibrate_arcframe.py

        nslits = len(rinput.slits_catalog)
        coeff_table = numpy.zeros((nslits, rinput.polynomial_degree + 1))

        image2d = hdulist[0].data
        naxis2, naxis1 = image2d.shape
        # read master table (TBM) and generate auxiliary parameters (valid for
        # all the slits) for the wavelength calibration
        wv_master = rinput.lines_catalog[:,0]
        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
          gen_triplets_master(wv_master, LDEBUG=True)
        # loop in number of slitlets
        for idx, slitdum in enumerate(rinput.slits_catalog):
            times_sigma = 3.0 # for minimum threshold level 
            nwinwidth=5 #number of pixels to detect and refine peaks (channels)
            # extract central spectrum
            xchannel, spectrum1d = slitdum.extract_midsp(image2d, nwidth=31,
                                                         method = 'mean',
                                                         LDEBUG = True)
            # find peaks (initial search providing integer numbers)
            threshold = numpy.median(spectrum1d)+times_sigma*sigmaG(spectrum1d)
            ipeaks_int = findPeaks_spectrum(spectrum1d, nwinwidth = nwinwidth, 
                                            data_threshold = threshold)
    
            # refine peaks fitting an appropriate function (providing float 
            # numbers)
            ipeaks_float = refinePeaks_spectrum(spectrum1d, ipeaks_int, nwinwidth, 
                                                method=1)
    
            # define interpolation function and interpolate the refined peak 
            # location, passing from index number (within the spectrum1d array) 
            # to channel number (note that this step takes care of the fact 
            # that the extracted spectrum may correspond to a subregion in the 
            # spectral direction)
            finterp_channel = interp1d(range(xchannel.size), xchannel, 
                                       kind='linear')
            xpeaks_refined = finterp_channel(ipeaks_float)

            # wavelength calibration
            solution = arccalibration_direct(wv_master,
                                             ntriplets_master,
                                             ratios_master_sorted,
                                             triplets_master_sorted_list,
                                             xpeaks_refined,
                                             naxis1,
                                             wv_ini_search = 8500, 
                                             wv_end_search = 25000,
                                             error_xpos_arc = 0.3,
                                             times_sigma_r = 3.0,
                                             frac_triplets_for_sum = 0.50,
                                             times_sigma_TheilSen = 10.0,
                                             poly_degree = 2,
                                             times_sigma_polfilt = 10.0,
                                             times_sigma_inclusion = 5.0,
                                             LDEBUG = False,
                                             LPLOT = False)

            numpy_array_with_coeff, crval1_approx, cdelt1_approx = \
              fit_solution(wv_master,
                           xpeaks_refined,
                           solution,
                           naxis1,
                           poly_degree = 2,
                           weighted = False,
                           LDEBUG = True,
                           LPLOT = True)

            coeff_table[idx] = numpy_array_with_coeff

        result = self.create_result(polynomial_coeffs=coeff_table)

        return result