#
# Copyright 2016-2017 Universidad Complutense de Madrid
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

"""
Spectroscopy mode, ABBA
"""

import astropy.io.fits as fits
import numina.core
import numina.exceptions
import numpy
from numina.array import combine
from numina.core import Product, RecipeError
from numina.core.requirements import ObservationResultRequirement

import emirdrp.datamodel
import emirdrp.decorators
import emirdrp.products as prods
from emirdrp.core import EmirRecipe
from emirdrp.processing.combine import basic_processing


class BaseABBARecipe(EmirRecipe):
    """Process images in ABBA mode"""

    obresult = ObservationResultRequirement()
    spec_abba = Product(prods.DataFrameType)

    # Accumulate 'spec_abba' results
    accum = Product(prods.DataFrameType, optional=True)

    @classmethod
    def build_recipe_input(cls, obsres, dal, pipeline='default'):
        return cls.build_recipe_input_gtc(obsres, dal, pipeline=pipeline)

    @classmethod
    def build_recipe_input_gtc(cls, obsres, dal, pipeline='default'):
        cls.logger.debug('start recipe input builder')
        stareImagesIds = obsres.stareSpectraIds
        cls.logger.debug('Stare Spectra images IDS: %s', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['stare'])

        naccum = obsres.naccum
        cls.logger.info('naccum: %d', naccum)
        if naccum != 1:  # if it is not the first dithering loop
            cls.logger.info("SEARCHING LATEST RESULT LS_ABBA TO ACCUMULATE")
            latest_result = dal.getLastRecipeResult("EMIR", "EMIR", "LS_ABBA")
            accum_dither = latest_result['elements']['accum']
            cls.logger.info("FOUND")
        else:
            cls.logger.info("NO ACCUMULATION LS_ABBA")
            accum_dither = stareImages[0]

        newOR = numina.core.ObservationResult()
        newOR.frames = stareImages
        newOR.naccum = naccum
        newOR.accum = accum_dither
        newRI = cls.create_input(obresult=newOR)
        cls.logger.debug('end recipe input builder')
        return newRI

    #@emirdrp.decorators.aggregate
    @emirdrp.decorators.loginfo
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
        result = self.create_result(spec_abba=hdulist)
        self.logger.info('end spectroscopy ABBA reduction')
        return result

    def process_abba(self, images):
        # Process for images in ABBA mode
        dataA0 = images[0][0].data
        dataB0 = images[1][0].data

        dataB1 = images[2][0].data
        dataA1 = images[3][0].data

        dataAB0 = dataA0 - dataB0
        dataAB1 = dataA1 - dataB1

        dataABBA = dataAB0 + dataAB1

        hdulist = self.create_proc_hdulist(images, dataABBA)
        self.logger.debug('update result header')
        hdu = hdulist[0]
        hdu.header['history'] = "Processed ABBA"
        hdu.header['NUM-NCOM'] = (2, 'Number of combined frames')
        dm = emirdrp.datamodel.EmirDataModel()
        for img, key in zip(images, ['A', 'B', 'B', 'A']):
            imgid = dm.get_imgid(img)
            hdu.header['history'] = "Image '{}' is '{}'".format(imgid, key)

        return hdulist

    def create_proc_hdulist(self, cdata, data_array):
        import astropy.io.fits as fits
        import uuid
        # Copy header of first image
        base_header = cdata[0][0].header.copy()

        hdu = fits.PrimaryHDU(data_array, header=base_header)
        self.set_base_headers(hdu.header)
        hdu.header['EMIRUUID'] = str(uuid.uuid1())
        # Update obsmode in header
        hdu.header['OBSMODE'] = 'LS_ABBA'
        # Headers of last image
        hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
        result = fits.HDUList([hdu])
        return result

    def create_accum_hdulist(self, cdata, data_array_n,
                             method_name='unkwnow', use_errors=False):
        import uuid

        base_header = cdata[0][0].header.copy()
        hdu = fits.PrimaryHDU(data_array_n[0], header=base_header)
        hdr = hdu.header
        self.set_base_headers(hdr)
        hdu.header['EMIRUUID'] = str(uuid.uuid1())
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
            num = fits.ImageHDU(data_array_n[2], name='MAP')
            hdulist = fits.HDUList([hdu, varhdu, num])
        else:
            hdulist = fits.HDUList([hdu])

        return hdulist

    def aggregate_result(self, partial_result, rinput):
        obresult = rinput.obresult
        # Check if this is our first run
        naccum = getattr(obresult, 'naccum', 0)
        accum = getattr(obresult, 'accum', None)
        # result to accumulate
        result_key = 'spec_abba'
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
