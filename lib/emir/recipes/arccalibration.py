#
# Copyright 2011-2015 Universidad Complutense de Madrid
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

from __future__ import division
from __future__ import print_function

'''Auxiliary Recipes for EMIR'''

import logging

import numpy

#from numina.core import BaseRecipeAutoQC
from numina.core import RecipeError
from numina.core import DataFrame#from numina.core import BaseRecipeAutoQC
from numina.core import Requirement, Product, Parameter
from numina.logger import log_to_history
from numina.array.combine import median
# from numina.flow.processing import BadPixelCorrector
from numina.core.requirements import ObservationResultRequirement
from numina.core.requirements import InstrumentConfigurationRequirement

import emir.instrument.channels as allchannels
from emir.core import EMIR_BIAS_MODES
from emir.core import gather_info_frames
from emir.core import EmirRecipe
from emir.dataproducts import MasterBias, MasterDark
from emir.dataproducts import MasterIntensityFlat
from emir.dataproducts import WavelengthCalibration, MasterSpectralFlat
from emir.dataproducts import ChannelLevelStatisticsType
from emir.dataproducts import ChannelLevelStatistics
from emir.dataproducts import SlitTransmissionCalibration
from emir.dataproducts import CoordinateList1DType
from emir.dataproducts import ArrayType
from emir.requirements import MasterBiasRequirement
from emir.requirements import MasterDarkRequirement
from emir.requirements import MasterBadPixelMaskRequirement
from emir.requirements import MasterIntensityFlatFieldRequirement
from emir.requirements import MasterSpectralFlatFieldRequirement
from .aiv.flows import init_filters_bdf
from .aiv.flows import init_filters_bd
from .aiv.flows import init_filters_b
from .aiv.flows import basic_processing_with_combination

from emir.wavecal.slitlet import Slitlet
from emir.wavecal.zscale import zscale
from emir.wavecal.statsummary import sigmaG
from emir.wavecal.arccalibration import arccalibration_direct, fit_solution, \
                                        gen_triplets_master
from emir.wavecal.findpeaks1D import findPeaks_spectrum, refinePeaks_spectrum
from scipy.interpolate import interp1d

#------------------------------------------------------------------------------

import numpy as np
from numpy.polynomial import polynomial
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------

_logger = logging.getLogger('numina.recipes.emir')


class ArcCalibrationRecipe(EmirRecipe):
    '''Bla, bla, bla.


    Bla bla bla
    '''

    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    polynomial_degree = Parameter(2, 'Polynomial degree of the arc calibration')

    polynomial_coeffs = Product(ArrayType)

    def run(self, rinput):
        _logger.info('starting arc calibration')

        flow = init_filters_bd(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)


        # ToDo
        # - replace TBM (To Be Modified) by proper code
        # - eliminate TBR code (To Be Removed)
        # - Set LDEBUG and LPLOT to False
        # - read slitletdef.txt to read slitlet definition for current image
        # - loop in nslits to excecute the code in calibrate_arcframe.py

        image2d = hdulist[0].data
        naxis2, naxis1 = image2d.shape
        # read master table (TBM) and generate auxiliary parameters (valid for
        # all the slits) for the wavelength calibration
        master_table = np.genfromtxt('../master_list.txt')
        wv_master = master_table[:,0]
        ntriplets_master, ratios_master_sorted, triplets_master_sorted_list = \
          gen_triplets_master(wv_master, LDEBUG = True) 
        # loop in number of slitlets
        nslitlets = 1
        for islit in range(nslitlets):
            # define slitlet
            slitdum = Slitlet(100,2000,11,100)
            times_sigma = 3.0 # for minimum threshold level 
            nwinwidth=5 #number of pixels to detect and refine peaks (channels)
            # extract central spectrum
            xchannel, spectrum1d = slitdum.extract_midsp(image2d, nwidth=31,
                                                         method = 'mean',
                                                         LDEBUG = True)
            # find peaks (initial search providing integer numbers)
            threshold = np.median(spectrum1d)+times_sigma*sigmaG(spectrum1d)
            ipeaks_int = findPeaks_spectrum(spectrum1d, nwinwidth = nwinwidth, 
                                            data_threshold = threshold,
                                            LDEBUG = True, LPLOT = False)
    
            # refine peaks fitting an appropriate function (providing float 
            # numbers)
            ipeaks_float = refinePeaks_spectrum(spectrum1d, ipeaks_int, nwinwidth, 
                                                method = 2, LDEBUG = False)
    
            # define interpolation function and interpolate the refined peak 
            # location, passing from index number (within the spectrum1d array) 
            # to channel number (note that this step takes care of the fact 
            # that the extracted spectrum may correspond to a subregion in the 
            # spectral direction)
            finterp_channel = interp1d(range(xchannel.size), xchannel, 
                                       kind = 'linear')
            xpeaks_refined = finterp_channel(ipeaks_float)
    
            if True: # TBR (to be removed in the future)
                print('xpeaks_refined (channel):\n',xpeaks_refined)
                # plot extracted spectrum
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_xlim([1,naxis1])
                ax.plot(xchannel,spectrum1d,'k-')
                ax.plot(xchannel[ipeaks_int], spectrum1d[ipeaks_int], 'ro')
                ax.plot([1,naxis1],[threshold, threshold], linestyle='dashed')
                ylim = ax.get_ylim()
                for xdum in zip(xpeaks_refined,xpeaks_refined):
                    ax.plot(xdum, ylim, linestyle='dotted', color='magenta')
                plt.show(block=False)
                raw_input('press <RETURN> to continue...')

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

            print('>>> approximate crval1, cdelt1:',crval1_approx,cdelt1_approx)
            print('>>> fitted coefficients.......:\n',numpy_array_with_coeff)
            raw_input('press <RETURN> to continue...')

        result = self.create_result(polynomial_coeffs = numpy_array_with_coeff)

        return result
