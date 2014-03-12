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
from numina.flow.processing import FlatFieldCorrector
from numina.flow import SerialFlow

from emir.core import RecipeResult
from emir.dataproducts import MasterBias, MasterDark, MasterBadPixelMask
from emir.dataproducts import FrameDataProduct, MasterIntensityFlat
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

        hdr.update('IMGTYP', 'BIAS', 'Image type')
        hdr.update('NUMTYP', 'MASTER_BIAS', 'Data product type')
        hdr.update('NUMXVER', __version__, 'Numina package version')
        hdr.update('NUMRNAM', self.__class__.__name__, 'Numina recipe name')
        hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

        exhdr = fits.Header()
        exhdr.update('extver', 1)
        varhdu = fits.ImageHDU(data[1], name='VARIANCE', header=exhdr)
        num = fits.ImageHDU(data[2], name='MAP')

        hdulist = fits.HDUList([hdu, varhdu, num])

        _logger.info('simple bias reduction ended')
 
        # qc is QC.UNKNOWN
        result = SimpleBiasRecipeResult(biasframe=DataFrame(hdulist))
        return result

def gather_info(hdulist):
    n_ext = len(hdulist)

    # READMODE is NUMERIC
    readmode = hdulist[0].header.get('READMODE', -1)
    readmods = hdulist[0].header.get('READMODS', 'undefined')
    bunit = hdulist[0].header.get('BUNIT', 'undefined')
    texp = hdulist[0].header.get('EXPTIME')
    adu_s = True
    if bunit:
        if bunit.lower() == 'adu':
            adu_s = False
        elif bunit.lower() == 'adu/s':
            adu_s = True
        else:
            _logger.warning('Unrecognized value for BUNIT %s', bunit)

    return {'n_ext': n_ext, 
            'readmode': readmode, 
            'readmods': readmods, 
            'adu_s': adu_s}


class TestBiasCorrectRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)

class TestBiasCorrectRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)

@define_requirements(TestBiasCorrectRecipeRequirements)
@define_result(TestBiasCorrectRecipeResult)
class TestBiasCorrectRecipe(BaseRecipe):

    def __init__(self):
        super(TestBiasCorrectRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple bias reduction')

        iinfo = []
        for frame in rinput.obresult.frames:
            with frame.open() as hdulist:
                iinfo.append(gather_info(hdulist))

        bias_info = {}

        if rinput.master_bias:
            with rinput.master_bias.open() as hdul:
                bias_info = gather_info(hdul)

        print(iinfo)
        print(bias_info)

        # SINGLE 0
        # CDS 1
        # FOWLER 2
        # RAMP 3
        # HDR_noseque $
        # BIAS 5

        for idx, ii in enumerate(iinfo):
            if not ii['readmode'].lower() in [0, 5]:
                # We have images in mode other than simple or bias BAD
                raise RecipeError('Image %d in inputs has READMODE %s' % (idx, ii['readmode']))
            #if not ii['readmode'].lower() in ['single', 'simple', 'bias']:
            #    # We have images in mode other than simple or bias BAD
            #    raise RecipeError('Image %d in inputs has READMODE %s', idx, ii.readmode)

        # Loading calibrations
        has_bias = False
        bias_corrector = None
        if rinput.master_bias:
            _logger.info('loading bias')
            has_bias = True
            with rinput.master_bias.open() as hdul:
                mbias = hdul[0].data.copy()
                bias_corrector = BiasCorrector(mbias)
        else:
            raise RecipeError("Bias required but not available")

        _logger.info('stacking %d images using median', len(cdata))
       
        cdata = []
        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open() # Check if I can return the same HDUList
                hdulist = bias_corrector(hdulist)
                cdata.append(hdulist)

            data = median([d['primary'].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)

        finally:
            for hdulist in cdata:
                hdulist.close()
            
        # Setup final header
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdulist = fits.HDUList([hdu])

        result = TestBiasCorrectRecipeResult(frame=hdulist)
        return result

class TestDarkCorrectRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')

class TestDarkCorrectRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)

@define_requirements(TestDarkCorrectRecipeRequirements)
@define_result(TestDarkCorrectRecipeResult)
class TestDarkCorrectRecipe(BaseRecipe):

    def __init__(self):
        super(TestDarkCorrectRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple dark reduction')

        iinfo = []
        for frame in rinput.obresult.frames:
            with frame.open() as hdulist:
                iinfo.append(gather_info(hdulist))

        bias_info = {}
        dark_info = {}

        if rinput.master_bias:
            with rinput.master_bias.open() as hdul:
                bias_info = gather_info(hdul)

        with rinput.master_dark.open() as hdul:
            dark_info = gather_info(hdul)

        print('images:', iinfo)
        print('bias:', bias_info)
        print('dark:', dark_info)

        # Loading calibrations
        if rinput.master_bias:
            _logger.info('loading bias')
            with rinput.master_bias.open() as hdul:
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbias)
        else:
            bias_corrector = IdNode()
            
        with rinput.master_dark.open() as mdark_hdul:
            _logger.info('loading dark')
            mdark = mdark_hdul[0].data
            dark_corrector = DarkCorrector(mdark)

        flow = SerialFlow([bias_corrector, dark_corrector])

        cdata = []
        try:
            for frame in rinput.obresult.frames:
                hdulist = frame.open() # Check if I can return the same HDUList
                hdulist = flow(hdulist)
                cdata.append(hdulist)

            data = median([d['primary'].data for d in cdata], dtype='float32')
            hdu = fits.PrimaryHDU(data[0], header=cdata[0]['primary'].header)

        finally:
            for hdulist in cdata:
                hdulist.close()
            
        hdr = hdu.header
        hdr['NUMXVER'] = (__version__, 'Numina package version')
        hdr['NUMRNAM'] = (self.__class__.__name__, 'Numina recipe name')
        hdulist = fits.HDUList([hdu])
        result = TestDarkCorrectRecipeResult(frame=hdulist)
        return result

class TestFlatCorrectRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')

class TestFlatCorrectRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)

@define_requirements(TestFlatCorrectRecipeRequirements)
@define_result(TestFlatCorrectRecipeResult)
class TestFlatCorrectRecipe(BaseRecipe):

    def __init__(self):
        super(TestFlatCorrectRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple flat reduction')

        # Loading calibrations
        if rinput.master_bias:
            _logger.info('loading bias')
            with rinput.master_bias.open() as hdul:
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbias)
        else:
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
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')
    master_sky = DataProductRequirement(MasterIntensityFlat, 'Master Sky calibration')

class TestSkyCorrectRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)

@define_requirements(TestSkyCorrectRecipeRequirements)
@define_result(TestSkyCorrectRecipeResult)
class TestSkyCorrectRecipe(BaseRecipe):

    def __init__(self):
        super(TestSkyCorrectRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple sky reduction')

        # Loading calibrations
        if rinput.master_bias:
            _logger.info('loading bias')
            with rinput.master_bias.open() as hdul:
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbias)
        else:
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


class TestPinholeRecipeRequirements(RecipeRequirements):
    obresult = ObservationResultRequirement()
    master_bias = DataProductRequirement(MasterBias, 'Master bias calibration', optional=True)
    master_dark = DataProductRequirement(MasterDark, 'Master dark calibration')
    master_flat = DataProductRequirement(MasterIntensityFlat, 'Master intensity flat calibration')
    master_sky = DataProductRequirement(MasterIntensityFlat, 'Master Sky calibration')
    pinhole_nominal_positions = None

class TestPinholetRecipeResult(RecipeResult):
    frame = Product(FrameDataProduct)

@define_requirements(TestPinholeRecipeRequirements)
@define_result(TestPinholeRecipeResult)
class TestPinholeRecipe(BaseRecipe):

    def __init__(self):
        super(TestPinholeRecipe, self).__init__(author=_s_author, 
            version="0.1.0")

    def run(self, rinput):
        _logger.info('starting simple sky reduction')

        # Loading calibrations
        if rinput.master_bias:
            _logger.info('loading bias')
            with rinput.master_bias.open() as hdul:
                mbias = hdul[0].data
                bias_corrector = BiasCorrector(mbias)
        else:
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


