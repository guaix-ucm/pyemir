#
# Copyright 2019 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Spectroscopy mode, combine AB or ABBA observations
"""

import astropy.io.fits as fits
import contextlib
import uuid

from numina.core import Result
from numina.core import Parameter
import numina.array.combine as combine

import emirdrp.requirements as reqs
from emirdrp.core.recipe import EmirRecipe
import emirdrp.datamodel
import emirdrp.products as prods
from numina.processing.combine import combine_imgs
from emirdrp.processing.wavecal.apply_rectwv_coeff import apply_rectwv_coeff
from emirdrp.processing.wavecal.median_slitlets_rectified import \
    median_slitlets_rectified
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import save_four_ds9
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import \
    save_spectral_lines_ds9


class ABBASpectraRectwv(EmirRecipe):
    """Process images in AB or ABBA  mode applying wavelength calibration

    Note that in this case the wavelength calibration has already been
    determined.

    """

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    rectwv_coeff = reqs.RectWaveCoeffRequirement()

    pattern = Parameter(
        'ABBA',
        description='Observation pattern',
        choices=['A', 'AB', 'ABBA']
    )
    repeat = Parameter(
        1,
        description='Repetitions at each individual A or B location'
    )
    nsequences = Parameter(
        1,
        description='Number of pattern sequences'
    )
    # do not allow 'sum' as combination method in order to avoid
    # confusion with number of counts in resulting image
    method = Parameter(
        'mean',
        description='Combination method',
        choices=['mean', 'median', 'sigmaclip']
    )

    reduced_mos_abba = Result(prods.ProcessedMOS)

    def run(self, rinput):
        self.logger.info('applying existing rect.+wavecal. calibration of '
                         'ABBA spectra')

        pattern = rinput.pattern
        repeat = rinput.repeat
        nsequences = rinput.nsequences

        self.logger.info(rinput.rectwv_coeff)
        self.logger.info('observation pattern......: {}'.format(pattern))
        self.logger.info('repeat...................: {}'.format(repeat))
        self.logger.info('nsequences...............: {}'.format(nsequences))

        # check number of images
        nimages = len(rinput.obresult.frames)
        pattern_sequence = ''
        for i in range(len(rinput.pattern)):
            pattern_sequence += pattern[i] * repeat
        self.logger.info('expected sequence pattern: {}'.format(
            pattern_sequence))
        if nimages % len(pattern_sequence) != 0:
            raise ValueError('Unexpected number of images in observation '
                             'result file')
        nexpected_images = len(pattern_sequence) * repeat * nsequences
        if nimages != nexpected_images:
            raise ValueError('Unexpected number of images in observation '
                             'result file')
        full_set = pattern_sequence * nsequences
        self.logger.info('full set of images.......: {}'.format(full_set))

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # available combination methods
        fmethod = getattr(combine, rinput.method)

        # basic reduction of A images
        list_a = [rinput.obresult.frames[i] for i, char in enumerate(full_set)
                  if char == 'A']
        with contextlib.ExitStack() as stack:
            self.logger.info('starting basic reduction of A images')
            hduls = [stack.enter_context(fname.open()) for fname in list_a]
            reduced_image_a = combine_imgs(
                hduls,
                method=fmethod,
                method_kwargs={},  # parameters for method
                errors=False,
                prolog=None
            )
        # ToDo: seguir aqui
        # - la idea es aplicar flow() individualmente a cada frame original
        # - hay que realizar rectificación y calibración en l.d.o.
        # - hay que enmascarar líneas de cielo y generar una imagen sin ellas
        # - obtener un perfil espacial (para rendijas predefinidas)
        # - hacer crosscorrelación para medir offsets con respecto a una
        # rendija(s) de referencia
        if False:
            print(type(hduls[0]))
            print(type(reduced_image_a))
            input("Stop here!")
        reduced_image_a = flow(reduced_image_a)
        hdr = reduced_image_a[0].header
        self.set_base_headers(hdr)
        self.save_intermediate_img(reduced_image_a, 'reduced_image_a.fits')

        if rinput.pattern == 'A':
            reduced_data = reduced_image_a[0].data.astype('float32')
        else:
            # basic reduction of B images
            list_b = [rinput.obresult.frames[i] for i, char in
                      enumerate(full_set) if char == 'B']
            with contextlib.ExitStack() as stack:
                self.logger.info('starting basic reduction of B images')
                hduls = [stack.enter_context(fname.open()) for fname in list_b]
                reduced_image_b = combine_imgs(
                    hduls,
                    method=fmethod,
                    method_kwargs={},  # parameters for method
                    errors=False,
                    prolog=None
                )
            reduced_image_b = flow(reduced_image_b)
            hdr = reduced_image_b[0].header
            self.set_base_headers(hdr)
            self.save_intermediate_img(reduced_image_b, 'reduced_image_b.fits')
            # computation of A-B
            data_a = reduced_image_a[0].data.astype('float32')
            data_b = reduced_image_b[0].data.astype('float32')
            reduced_data = data_a - data_b

        # create reduced_image
        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(fname.open()) for fname in
                     rinput.obresult.frames]
            reduced_image = self.create_reduced_image(
                hduls,
                reduced_data,
                rinput.pattern,
                full_set
            )

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # apply rectification and wavelength calibration
        reduced_mos_abba = apply_rectwv_coeff(
            reduced_image,
            rinput.rectwv_coeff
        )

        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(rinput.rectwv_coeff)
            save_spectral_lines_ds9(rinput.rectwv_coeff)

        # compute median spectra employing the useful region of the
        # rectified image
        if self.intermediate_results:
            for imode, outfile in enumerate(['median_spectra_full',
                                             'median_spectra_slitlets',
                                             'median_spectrum_slitlets']):
                median_image = median_slitlets_rectified(
                    reduced_mos_abba, mode=imode
                )
                self.save_intermediate_img(median_image, outfile + '.fits')

        # save results in results directory
        self.logger.info('end rect.+wavecal. reduction of ABBA spectra')
        result = self.create_result(reduced_mos_abba=reduced_mos_abba)
        return result

    def create_reduced_image(self, hduls, reduced_data, pattern, full_set):
        # Copy header of first image
        base_header = hduls[0][0].header.copy()

        hdu = fits.PrimaryHDU(reduced_data, header=base_header)
        self.set_base_headers(hdu.header)

        self.logger.debug('update result header')
        hdu.header['UUID'] = str(uuid.uuid1())
        hdu.header['OBSMODE'] = pattern + ' pattern'
        hdu.header['TSUTC2'] = hduls[-1][0].header['TSUTC2']
        hdu.header['history'] = "Processed " + pattern + " pattern"
        hdu.header['NUM-NCOM'] = (len(hduls), 'Number of combined frames')
        dm = emirdrp.datamodel.EmirDataModel()
        for img, key in zip(hduls, full_set):
            imgid = dm.get_imgid(img)
            hdu.header['history'] = "Image '{}' is '{}'".format(imgid, key)

        result = fits.HDUList([hdu])
        return result

    def set_base_headers(self, hdr):
        newhdr = super(ABBASpectraRectwv, self).set_base_headers(hdr)
        # Update EXP to 0
        newhdr['EXP'] = 0
        return newhdr
