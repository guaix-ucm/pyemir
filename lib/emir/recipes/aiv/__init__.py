#
# Copyright 2013 Universidad Complutense de Madrid
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

'''AIV Recipes for EMIR'''

import logging

import numpy
import pyfits
from scipy.stats import linregress
import matplotlib.pyplot as plt

from numina.core import RecipeError
from numina.core import BaseRecipe, RecipeRequirements, DataFrame
from numina.core import Requirement, Product, DataProductRequirement
from numina.core import define_requirements, define_result
from numina.logger import log_to_history
from numina.array.combine import median
from numina import __version__
from numina.flow.processing import BiasCorrector

from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import NonLinearityCalibration, MasterIntensityFlat
from emir.dataproducts import DarkCurrentValue

_logger = logging.getLogger('numina.recipes.emir')

_s_author = "Sergio Pascual <sergiopr@fis.ucm.es>"
            
class DarkCurrentRecipeRequirements(RecipeRequirements):
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)

class DarkCurrentRecipeResult(RecipeResult):
    darkcurrent = Product(DarkCurrentValue)

@define_requirements(DarkCurrentRecipeRequirements)
@define_result(DarkCurrentRecipeResult)
class DarkCurrentRecipe(BaseRecipe):
    '''Recipe to process data taken in Dark Current image Mode.''' 

    taskidx = '2.1.03'
    taskname = 'Dark current (+contribution from instrument)'

    def __init__(self):
        super(DarkCurrentRecipe, self).__init__(author=_s_author, 
                    version="0.1.0")
        
    def run(self, obresult, reqs):
        _logger.info('starting dark current reduction')

        if reqs.master_bias is not None:
            _logger.debug("master bias=%s",reqs.master_bias)
            master_bias = pyfits.getdata(reqs.master_bias.filename)
            master_bias_base = master_bias
        else:
            master_bias_base = 0
            
        values_t = []
        values_d = []
        for frame in obresult.frames:
            with pyfits.open(frame.label) as hdulist:
                # FIXME: single images must be corrected to have an uniform 
                # exposure time
                texp = hdulist[0].header['exposed']
                corrected = hdulist[0].data - master_bias_base
                corrected_mean = corrected.mean()
                _logger.debug("%s mean=%f exposure=%8.2f", frame.label, corrected_mean, texp)
                values_t.append(texp)
                values_d.append(corrected_mean)

        values_t = numpy.array(values_t)
        values_d = numpy.array(values_d)
        slope, intercept, r_value, p_value, std_err = linregress(values_t, values_d)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Exposure time')
        ax.set_ylabel('Dark current [ADU]')
        ax.plot(values_t, values_d, '-*')
        ax.plot(values_t, slope * values_t + intercept, 'r-')
        fig.savefig('dark-current.png')
        print('slope=', slope, 'intercept=', intercept)

        _logger.info('dark current reduction ended')
        result = DarkCurrentRecipeResult(darkcurrent=DarkCurrentValue())
        return result

class SimpleBiasRecipeRequirements(RecipeRequirements):
    pass

class SimpleBiasRecipeResult(RecipeResult):
    biasframe = Product(MasterBias)

@define_requirements(SimpleBiasRecipeRequirements)
@define_result(SimpleBiasRecipeResult)
class SimpleBiasRecipe(BaseRecipe):
    '''    
    Recipe to process data taken in SimpleBias image Mode.

    Bias images only appear in Simple Readout mode.

    **Outputs:**

     * A combined bias frame, with variance extension.

    **Procedure:**
    
    The list of images can be readly processed by combining them 
    with a median algorithm.
    '''

    def __init__(self):
        super(SimpleBiasRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, obresult, reqs):
        _logger.info('starting simple bias reduction')

        cdata = []

        try:
            for frame in obresult.frames:
                hdulist = pyfits.open(frame.label, memmap=True, mode='readonly')
                cdata.append(hdulist)

            _logger.info('stacking %d images using median', len(cdata))
            
            data = median([d['primary'].data for d in cdata], dtype='float32')

            hdu = pyfits.PrimaryHDU(data[0], header=cdata[0]['PRIMARY'].header)

        finally:
            for hdulist in cdata:
                hdulist.close()
            
        # update hdu header with
        # reduction keywords
        hdr = hdu.header

        hdr.update('IMGTYP', 'BIAS', 'Image type')
        hdr.update('NUMTYP', 'MASTER_BIAS', 'Data product type')
        hdr.update('NUMXVER', __version__, 'Numina package version')
        hdr.update('NUMRNAM', self.__class__.__name__, 'Numina recipe name')
        hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

        exhdr = pyfits.Header()
        exhdr.update('extver', 1)
        varhdu = pyfits.ImageHDU(data[1], name='VARIANCE', header=exhdr)
        num = pyfits.ImageHDU(data[2], name='MAP')

        hdulist = pyfits.HDUList([hdu, varhdu, num])

        _logger.info('simple bias reduction ended')
 
        result = SimpleBiasRecipeResult(biasframe=DataFrame(hdulist))
        return result

