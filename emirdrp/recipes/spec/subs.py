#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Spectroscopy mode, ABBA
"""
import uuid
import pkgutil

import astropy.io.fits as fits
import numpy

import numina.core
import numina.exceptions
import numina.core.query as qmod
import numina.ext.gtc
from numina.array import combine
from numina.core import Result, Requirement
from numina.exceptions import RecipeError
from numina.core.requirements import ObservationResultRequirement
from numina.frame.utils import copy_img

import emirdrp.datamodel
import emirdrp.decorators
import emirdrp.products as prods
from emirdrp.core.recipe import EmirRecipe
from emirdrp.processing.combine import basic_processing


class BaseABBARecipe(EmirRecipe):
    """Process images in ABBA mode"""

    obresult = ObservationResultRequirement(
        query_opts=qmod.ResultOf(
            'STARE_SPECTRA.reduced_mos',
            node='children',
            id_field="stareSpectraIds"
        )
    )
    accum_in = Requirement(
        prods.DataFrameType,
        description='Accumulated result',
        optional=True,
        destination='accum',
        query_opts=qmod.ResultOf(
            'LS_ABBA.accum',
            node='prev'
        )
    )

    reduced_mos_abba = Result(prods.ProcessedMOS)
    # Accumulate 'reduced_mos_abba' results
    accum = Result(prods.ProcessedMOS, optional=True)

    def build_recipe_input(self, obsres, dal, pipeline='default'):
        if numina.ext.gtc.check_gtc():
            self.logger.debug('running in GTC environment')
            return self.build_recipe_input_gtc(obsres, dal)
        else:
            self.logger.debug('running outside of GTC environment')
            return super(BaseABBARecipe, self).build_recipe_input(
                obsres, dal
            )

    def build_recipe_input_gtc(self, obsres, dal):
        self.logger.debug('start recipe input builder')
        stareImagesIds = obsres.stareSpectraIds
        self.logger.debug('Stare Spectra images IDS: %s', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['reduced_mos'])

        naccum = obsres.naccum
        self.logger.info('naccum: %d', naccum)
        if naccum != 1:  # if it is not the first dithering loop
            self.logger.info("SEARCHING LATEST RESULT LS_ABBA TO ACCUMULATE")
            latest_result = dal.getLastRecipeResult("EMIR", "EMIR", "LS_ABBA")
            accum_dither = latest_result['elements']['accum']
            self.logger.info("FOUND")
        else:
            self.logger.info("NO ACCUMULATION LS_ABBA")
            accum_dither = stareImages[0]

        newOR = numina.core.ObservationResult()
        newOR.frames = stareImages
        newOR.naccum = naccum
        newOR.accum = accum_dither
        newRI = self.create_input(obresult=newOR)
        self.logger.debug('end recipe input builder')
        return newRI

    def run(self, rinput):
        partial_result = self.run_single(rinput)
        new_result = self.aggregate_result(partial_result, rinput)
        return new_result

    def run_single(self, rinput):
        self.logger.info('starting spectroscopy ABBA reduction')

        flow = self.init_filters(rinput)
        nimages = len(rinput.obresult.frames)
        self.logger.info('we receive %d images', nimages)
        if nimages != 4:
            msg = 'Recipe expects 4 images, received %d' % nimages
            raise numina.exceptions.RecipeError(msg)

        procesed_hdulists = basic_processing(rinput, flow)

        # INPUTS are ABBA, so
        #
        hdulist = self.process_abba(procesed_hdulists)
        grism = hdulist[0].header.get('GRISM', 'unknown')
        if grism.lower() == 'open':
            # perform TEST10 in addition
            import emirdrp.recipes.acquisition.maskcheck as mk
            from numina.core import ObservationResult, DataFrame
            import numpy
            try:
                import StringIO as S
            except ImportError:
                import io as S

            self.logger.info('GRISM is OPEN, doing a RECIPE10')
            sub = mk.MaskCheckRecipe()
            sub.configure(instrument='EMIR', mode='TEST10')
            o = ObservationResult()
            o.__dict__ = rinput.obresult.__dict__
            o.frames = [DataFrame(frame=hdulist)]
            subd = {}
            subd['obresult'] = o
            bardata = pkgutil.get_data('emirdrp.instrument.configs', 'bars_nominal_positions_test.txt')
            ss = S.StringIO(bardata.decode('utf8'))
            subd['bars_nominal_positions'] = numpy.loadtxt(ss)

            subinput = mk.MaskCheckRecipe.RecipeInput(**subd)

            sub.run(subinput)

        result = self.create_result(reduced_mos_abba=hdulist)
        self.logger.info('end spectroscopy ABBA reduction')
        return result

    def process_abba(self, images):
        """Process four images in ABBA mode"""
        # FIXME: Duplicated in processing
        dataA0 = images[0][0].data.astype('float32')
        dataB0 = images[1][0].data.astype('float32')

        dataB1 = images[2][0].data.astype('float32')
        dataA1 = images[3][0].data.astype('float32')

        dataAB0 = dataA0 - dataB0
        dataAB1 = dataA1 - dataB1

        dataABBA = dataAB0 + dataAB1

        hdulist = self.create_proc_hdulist(images, dataABBA)
        self.logger.debug('update result header')
        hdu = hdulist[0]
        hdu.header['history'] = "Processed ABBA"
        hdu.header['NUM-NCOM'] = (2, 'Number of combined frames')
        hdu.header['IMGOBBL'] = 0
        dm = emirdrp.datamodel.EmirDataModel()
        for img, key in zip(images, ['A', 'B', 'B', 'A']):
            imgid = dm.get_imgid(img)
            hdu.header['history'] = "Image '{}' is '{}'".format(imgid, key)

        return hdulist

    def create_proc_hdulist(self, cdata, data_array):
        # Copy header of first image
        result = copy_img(cdata[0])

        hdu = result[0]
        hdu.data = data_array
        self.set_base_headers(hdu.header)
        hdu.header['UUID'] = str(uuid.uuid1())
        # Update obsmode in header
        hdu.header['OBSMODE'] = 'LS_ABBA'
        # Headers of last image
        hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
        return result

    def create_accum_hdulist(self, cdata, data_array_n,
                             method_name='unkwnow', use_errors=False):
        # FIXME: duplicated
        result = copy_img(cdata[0])
        hdu = result[0]
        hdu.data = data_array_n[0]
        hdr = hdu.header
        self.set_base_headers(hdr)
        hdu.header['UUID'] = str(uuid.uuid1())
        hdr['IMGOBBL'] = 0
        hdr['TSUTC2'] = cdata[-1][0].header['TSUTC2']

        hdu.header['history'] = "Combined %d images using '%s'" % (
            len(cdata),
            method_name
        )
        #hdu.header['history'] = 'Combination time {}'.format(
        #    datetime.datetime.utcnow().isoformat()
        #)
        # Update NUM-NCOM, sum of individual frames
        ncom = 0
        for hdul in cdata:
            ncom += hdul[0].header['NUM-NCOM']
        hdr['NUM-NCOM'] = ncom
        #
        if use_errors:
            varhdu = fits.ImageHDU(data_array_n[1], name='VARIANCE')
            result.append(varhdu)
            num = fits.ImageHDU(data_array_n[2], name='MAP')
            result.append(num)
        return result

    def aggregate_result(self, partial_result, rinput):
        obresult = rinput.obresult
        # Check if this is our first run
        naccum = getattr(obresult, 'naccum', 0)
        accum = getattr(obresult, 'accum', None)
        # result to accumulate
        result_key = 'reduced_mos_abba'
        field_to_accum = getattr(partial_result, result_key)

        if naccum == 0:
            self.logger.debug('naccum is not set, do not accumulate')
            return partial_result
        elif naccum == 1:
            self.logger.debug('round %d initialize accumulator', naccum)
            newaccum = field_to_accum
        elif naccum > 1:
            self.logger.debug('round %d of accumulation', naccum)
            newaccum = self.aggregate_frames(accum, field_to_accum, naccum)
        else:
            msg = 'naccum set to %d, invalid' % (naccum,)
            self.logger.error(msg)
            raise RecipeError(msg)

        # Update partial result
        partial_result.accum = newaccum

        return partial_result

    def aggregate_frames(self, accum, frame, naccum):
        return self.aggregate2(accum, frame, naccum)

    def aggregate2(self, img1, img2, naccum):

        frames = [img1, img2]
        use_errors = True
        # Initial checks
        fframe = frames[0]
        # Ref image
        img = fframe.open()
        has_num_ext = 'NUM' in img
        has_bpm_ext = 'BPM' in img
        base_header = img[0].header
        baseshape = img[0].shape

        data_hdul = []
        for f in frames:
            img = f.open()
            data_hdul.append(img)

        if has_num_ext:
            self.logger.debug('Using NUM extension')
            masks = [numpy.where(m['NUM'].data, 0, 1).astype('uint8') for m in data_hdul]
        elif has_bpm_ext:
            self.logger.debug('Using BPM extension')
            masks = [m['BPM'].data for m in data_hdul]
        else:
            self.logger.warning('BPM missing, use zeros instead')
            false_mask = numpy.zeros(baseshape, dtype='int16')
            masks = [false_mask for _ in data_hdul]

        self.logger.info('Combine target images (final, aggregate)')

        weight_accum = 2 * (1 - 1.0 / naccum)
        weight_frame = 2.0 / naccum
        self.logger.debug("weights for 'accum' and 'frame', %s", [weight_accum, weight_frame])
        scales = [1.0 / weight_accum, 1.0 / weight_frame]
        method = combine.mean
        data_arr = [hdul[0].data for hdul in data_hdul]
        out = method(data_arr, masks=masks, scales=scales, dtype='float32')

        self.logger.debug('create result image')

        return self.create_accum_hdulist(
            data_hdul,
            out,
            method_name=method.__name__,
            use_errors=False
        )
