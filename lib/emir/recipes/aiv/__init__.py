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

'''AIV Recipes for EMIR'''

import logging

import numpy
from astropy.io import fits
from scipy.stats import linregress
import matplotlib.pyplot as plt

from numina.core import RecipeError
from numina.core import BaseRecipe, RecipeRequirements, DataFrame
from numina.core import Requirement, Product, DataProductRequirement
from numina.core import define_requirements, define_result
from numina.core.requirements import ObservationResultRequirement

from numina.array.combine import median, mean
from numina import __version__
from numina.flow.processing import BiasCorrector, DarkCorrector
from numina.flow.processing import FlatFieldCorrector, SkyCorrector
from numina.flow.node import IdNode
from numina.flow import SerialFlow

from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import DataFrameType, MasterIntensityFlat
from emir.dataproducts import DarkCurrentValue, CoordinateList2DType
from emir.core import gather_info
from emir.core import EMIR_BIAS_MODES

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
            master_bias = fits.getdata(reqs.master_bias.filename)
            master_bias_base = master_bias
        else:
            master_bias_base = 0
            
        values_t = []
        values_d = []
        for frame in obresult.frames:
            with fits.open(frame.label) as hdulist:
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
        slope, intercept, _r_value, _p_value, _std_err = linregress(values_t, values_d)
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
    obresult = ObservationResultRequirement()

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

    def run(self, inputs):
        _logger.info('starting simple bias reduction')


        # On Pyhton >= 3.3 we can use ExitStack as context manager
        # for all the file openings
        cdata = []
        try:
            for frame in inputs.obresult.frames:
                cdata.append(frame.open())

            _logger.info('stacking %d images using median', len(cdata))
            
            data = median([d['primary'].data for d in cdata], dtype='float32')
            template_head = cdata[0]['PRIMARY'].header
            hdu = fits.PrimaryHDU(data[0], header=template_head)

        finally:
            for hdulist in cdata:
                hdulist.close()
            
        # update hdu header with
        # reduction keywords
        hdr = hdu.header

        hdr['IMGTYP'] = ('BIAS', 'Image type')
        hdr['NUMTYP'] = ('MASTER_BIAS', 'Data product type')
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')

        exhdr = fits.Header()
        exhdr['extver'] = 1
        varhdu = fits.ImageHDU(data[1], name='VARIANCE', header=exhdr)
        num = fits.ImageHDU(data[2], name='MAP')

        hdulist = fits.HDUList([hdu, varhdu, num])

        _logger.info('simple bias reduction ended')
 
        # qc is QC.UNKNOWN
        result = SimpleBiasRecipeResult(biasframe=DataFrame(hdulist))
        return result




class TestBiasCorrectRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration')

class TestBiasCorrectRecipeResult(RecipeResult):
    frame = Product(DataFrameType)

@define_requirements(TestBiasCorrectRecipeRequirements)
@define_result(TestBiasCorrectRecipeResult)
class TestBiasCorrectRecipe(BaseRecipe):

    def __init__(self):
        super(TestBiasCorrectRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple bias reduction')
        
        meta = gather_info(rinput)
        iinfo = meta['obresult']
        
        use_bias = False
        

        mode = iinfo[0]['readmode']
        if mode.lower() in EMIR_BIAS_MODES:
            _logger.info('readmode is %s, bias required', mode)
            use_bias = True
                
        else:
            _logger.error('readmode is %s, no bias required', mode)
            raise RecipeError('readmode is %s, no bias required', mode)
             
        bias_info = meta['master_bias']
           
        print('images info:', iinfo)
        print('bias info:', bias_info)

        # Loading calibrations
        if use_bias:
            with rinput.master_bias.open() as hdul:
                _logger.info('loading bias')
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbias)
        else:
            _logger.info('ignoring bias')
            bias_corrector = IdNode()

        flow = SerialFlow([bias_corrector])
       
        cdata = []
        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                hdulist = flow(hdulist)
                cdata.append(hdulist)

            _logger.info('combining images using median')
            data = median([d['primary'].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)

        finally:
            for hdulist in cdata:
                hdulist.close()
            
        _logger.info('stacking %d images using median', len(cdata))
        # Setup final header
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdulist = fits.HDUList([hdu])

        result = TestBiasCorrectRecipeResult(frame=hdulist)
        return result

class TestDarkCorrectRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration')
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')

class TestDarkCorrectRecipeResult(RecipeResult):
    frame = Product(DataFrameType)

@define_requirements(TestDarkCorrectRecipeRequirements)
@define_result(TestDarkCorrectRecipeResult)
class TestDarkCorrectRecipe(BaseRecipe):

    def __init__(self):
        super(TestDarkCorrectRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple dark reduction')

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
                
        bias_info = meta['master_bias']
        dark_info = meta['master_dark']

        print('images info:', iinfo)
        if use_bias:
            print('bias info:', bias_info)
        print('dark info:', dark_info)

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

        flow = SerialFlow([bias_corrector, dark_corrector])

        cdata = []
        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                hdulist = flow(hdulist)
                cdata.append(hdulist)

            _logger.info('stacking %d images using median', len(cdata))
            data = median([d['primary'].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)

        finally:
            for hdulist in cdata:
                hdulist.close()
            
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdulist = fits.HDUList([hdu])
        
        result = TestDarkCorrectRecipeResult(frame=hdulist)
        
        return result

class TestFlatCorrectRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration')
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')

class TestFlatCorrectRecipeResult(RecipeResult):
    frame = Product(DataFrameType)

@define_requirements(TestFlatCorrectRecipeRequirements)
@define_result(TestFlatCorrectRecipeResult)
class TestFlatCorrectRecipe(BaseRecipe):

    def __init__(self):
        super(TestFlatCorrectRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple flat reduction')

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
        
        bias_info = meta['master_bias']
        dark_info = meta['master_dark']
        flat_info = meta['master_flat']


        print('images info:', iinfo)
        if use_bias:
            print('bias info:', bias_info)
        print('dark info:', dark_info)
        print('flat info:', flat_info)

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

        flow = SerialFlow([bias_corrector, dark_corrector, flat_corrector])

        cdata = []
        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                final = flow(hdulist)
                cdata.append(final)

            _logger.info('stacking %d images using median', len(cdata))
            data = median([d['primary'].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)

        finally:
            for hdulist in cdata:
                hdulist.close()
            
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        hdulist = fits.HDUList([hdu])
        result = TestFlatCorrectRecipeResult(frame=hdulist)

        return result


class TestSkyCorrectRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration')
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')
    master_sky = DataProductRequirement(MasterIntensityFlat, 'Master Sky calibration')

class TestSkyCorrectRecipeResult(RecipeResult):
    frame = Product(DataFrameType)

@define_requirements(TestSkyCorrectRecipeRequirements)
@define_result(TestSkyCorrectRecipeResult)
class TestSkyCorrectRecipe(BaseRecipe):

    def __init__(self):
        super(TestSkyCorrectRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple sky reduction')

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

        cdata = []
        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open()
                final = flow(hdulist)
                cdata.append(final)

            _logger.info('stacking %d images using median', len(cdata))
            data = median([d['primary'].data for d in cdata], dtype='float32')
            _logger.debug('create result Primary HDU')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)

        finally:
            _logger.debug('closing images')
            for hdulist in cdata:
                hdulist.close()
            
        _logger.debug('update result header')
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdr['NUMRVER'] = (self.__version__, 'Numina recipe version')
        # Update SEC to 0
        hdr['SEC'] = 0
        hdulist = fits.HDUList([hdu])

        result = TestSkyCorrectRecipeResult(frame=hdulist)

        return result
