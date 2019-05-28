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

from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
import contextlib
import numpy as np
from scipy import ndimage
import uuid

from numina.array.wavecalib.crosscorrelation import periodic_corr1d
from numina.core import Result
from numina.core import Parameter
import numina.array.combine as combine
from numina.array.distortion import shift_image2d

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
from emirdrp.processing.wavecal.useful_mos_xpixels import useful_mos_xpixels

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2


def get_isky(i, pattern):
    """Return index number of image to be employed as sky"""

    if pattern == 'AB':
        if (i + 1) % 2 == 0:
            return i - 1
        else:
            return i + 1
    elif pattern == 'ABBA':
        n_previous_sequences = (i + 1) // 4
        if (i + 1) % 4 == 0:
            n_previous_sequences -= 1
        ieff = i - n_previous_sequences * 4
        if ieff == 0:
            iout = i + 1
        elif ieff == 1:
            iout = i - 1
        elif ieff == 2:
            iout = i + 1
        elif ieff == 3:
            iout = i - 1
        else:
            raise ValueError('Unexpected ieff value={}'.format(ieff))
    else:
        raise ValueError('Unexpected pattern: {}'.format(pattern))

    return iout


def compute_wcs_offsets(hduls):
    """Compute offsets between ABBA images from WCS info in image headers"""

    nimages = len(hduls)
    c0 = None
    list_sep_arcsec = []
    for i in range(nimages):
        wcsh = wcs.WCS(hduls[i][0].header)
        ra, dec = wcsh.wcs_pix2world(EMIR_NAXIS1 / 2 + 0.5,
                                     EMIR_NAXIS2 / 2 + 0.5, 1)
        if i == 0:
            c0 = SkyCoord(ra, dec, unit="deg")
        c = SkyCoord(ra, dec, unit="deg")
        sep = c.separation(c0).arcsec
        list_sep_arcsec.append(round(sep, 4))

    return list_sep_arcsec


class ABBASpectraRectwv(EmirRecipe):
    """Process images in AB or ABBA  mode applying wavelength calibration

    Note that in this case the wavelength calibration has already been
    determined.

    """

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    rectwv_coeff = reqs.RectWaveCoeffRequirement(optional=True)
    list_rectwv_coeff = reqs.ListOfRectWaveCoeffRequirement(optional=True)

    pattern = Parameter(
        'ABBA',
        description='Observation pattern',
        choices=['AB', 'ABBA']
    )
    # do not allow 'sum' as combination method in order to avoid
    # confusion with number of counts in resulting image
    method = Parameter(
        'mean',
        description='Combination method',
        choices=['mean', 'median', 'sigmaclip']
    )
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method',
        optional=True
    )
    refine_target_along_slitlet = Parameter(
        dict(),
        description='Parameters to refine location of target along the slitlet'
    )

    reduced_mos_abba = Result(prods.ProcessedMOS)

    def run(self, rinput):
        nimages = len(rinput.obresult.frames)
        pattern = rinput.pattern
        pattern_length = len(pattern)

        # check pattern sequence matches number of images
        if nimages % pattern_length != 0:
            raise ValueError('Number of images is not a multiple of pattern '
                             'length: {}, {}'.format(nimages, pattern_length))
        nsequences = nimages // pattern_length

        # check rectification and wavelength calibration information
        if rinput.rectwv_coeff is None and rinput.list_rectwv_coeff is None:
            raise ValueError('No rectwv_coeff nor list_rectwv_coeff data have '
                             'been provided')
        elif rinput.rectwv_coeff is not None and \
                rinput.list_rectwv_coeff is not None:
            raise ValueError("rectwv_coeff and list_rectwv_coeff cannot be "
                             "used simultaneously")
        elif rinput.rectwv_coeff is not None:
            list_rectwv_coeff = [rinput.rectwv_coeff] * nimages
        elif rinput.list_rectwv_coeff is not None:
            if len(rinput.list_rectwv_coeff) != nimages:
                raise ValueError("Unexpected number of rectwv_coeff files "
                                 "in list_rectwv_coeff")
            else:
                list_rectwv_coeff = rinput.list_rectwv_coeff
                # check filter and grism are the same in all JSON files
                for item in ['grism', 'filter']:
                    list_values = []
                    for calib in list_rectwv_coeff:
                        list_values.append(calib.tags[item])
                if len(set(list_values)) != 1:
                    raise ValueError('list_rectwv_coeff contains coefficients'
                                     'for different {}s')
        else:
            raise ValueError("Unexpected error!")

        grism_name = list_rectwv_coeff[0].tags['grism']
        filter_name = list_rectwv_coeff[0].tags['filter']

        # compute offsets from WCS info in image headers
        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(fname.open()) for fname in
                     rinput.obresult.frames]
            list_sep_arcsec = compute_wcs_offsets(hduls)

        for i in range(nimages):
            self.logger.info(list_rectwv_coeff[i])
        self.logger.info('observation pattern.......: {}'.format(pattern))
        self.logger.info('nsequences................: {}'.format(nsequences))
        self.logger.info('offsets (arcsec, from WCS): {}'.format(
            list_sep_arcsec))

        full_set = pattern * nsequences
        self.logger.info('full set of images.: {}'.format(full_set))

        # basic parameters to determine useful pixels in the wavelength
        # direction
        dumdict = rinput.refine_target_along_slitlet

        if 'vertical_pix_region_a' in dumdict.keys():
            vertical_pix_region_a = dumdict['vertical_pix_region_a']
        else:
            vertical_pix_region_a = None

        if 'vertical_pix_region_b' in dumdict.keys():
            vertical_pix_region_b = dumdict['vertical_pix_region_b']
        else:
            vertical_pix_region_b = None

        try:
            npix_removed_around_ohlines = \
                int(dumdict['npix_removed_around_ohlines'])
        except KeyError:
            npix_removed_around_ohlines = 0
        except ValueError:
            raise ValueError(
                'wrong value: npix_removed_around_ohlines={}'.format(
                dumdict['npix_removed_around_ohlines'])
            )

        if 'list_valid_wvregions_a' in dumdict.keys():
            if vertical_pix_region_a is None:
                raise ValueError('Unexpected list_valid_wvregions_a when '
                                 'vertical_pix_region_a is not set')
            list_valid_wvregions_a = dumdict['list_valid_wvregions_a']
        else:
            list_valid_wvregions_a = None

        if 'list_valid_wvregions_b' in dumdict.keys():
            if vertical_pix_region_b is None:
                raise ValueError('Unexpected list_valid_wvregions_b when '
                                 'vertical_pix_region_b is not set')

            list_valid_wvregions_b = dumdict['list_valid_wvregions_b']
        else:
            list_valid_wvregions_b = None

        try:
            nwidth_medfilt = int(dumdict['nwidth_medfilt'])
        except KeyError:
            nwidth_medfilt = 0
        except ValueError:
            raise ValueError('wrong value: nwidth_medfilt={}'.format(
                dumdict['nwidth_medfilt']))

        self.logger.info('vertical_pix_region_a: {}'.format(
            vertical_pix_region_a))
        self.logger.info('vertical_pix_region_b: {}'.format(
            vertical_pix_region_b))
        self.logger.info('npix_removed_around_ohlines: {}'.format(
            npix_removed_around_ohlines))
        self.logger.info('list_valid_wvregions_a: {}'.format(
            list_valid_wvregions_a))
        self.logger.info('list_valid_wvregions_b: {}'.format(
            list_valid_wvregions_b))
        self.logger.info('nwidth_medfilt: {}'.format(nwidth_medfilt))

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # basic reduction of A images
        self.logger.info('starting reduction of individual images')
        first_a = True
        first_b = True
        xisok_a = None
        xisok_b = None
        reference_profile_a = None
        reference_profile_b = None
        refine_a = vertical_pix_region_a is not None
        refine_b = vertical_pix_region_b is not None
        list_offsets = []
        for i, char in enumerate(full_set):
            if char == 'A' and refine_a:
                refine_image = True
            elif char == 'B' and refine_b:
                refine_image = True
            else:
                refine_image = False
            if refine_image:
                frame = rinput.obresult.frames[i]
                isky = get_isky(i, pattern)
                frame_sky = rinput.obresult.frames[isky]
                self.logger.info('image {} ({} of {})'.format(char, i + 1,
                                                              nimages))
                self.logger.info('image: {}'.format(frame.filename))
                self.logger.info('(sky): {}'.format(frame_sky.filename))
                with frame.open() as f:
                    base_header = f[0].header.copy()
                    # check grism
                    grism_name_ = base_header['grism']
                    if grism_name_ != grism_name:
                        raise ValueError('Unexpected grism: {}'.format(
                            grism_name_))
                    # check filter
                    filter_name_ = base_header['filter']
                    if filter_name_ != filter_name:
                        raise ValueError('Unexpected filter: {}'.format(
                            filter_name_))
                    data = f[0].data.astype('float32')
                with frame_sky.open() as fsky:
                    data_sky = fsky[0].data.astype('float32')
                data -= data_sky
                hdu = fits.PrimaryHDU(data, header=base_header)
                hdu.header['UUID'] = str(uuid.uuid1())
                hdul = fits.HDUList([hdu])
                for hdu in f[1:]:
                    hdul.append(hdu.copy())
                # basic reduction
                reduced_image = flow(hdul)
                hdr = reduced_image[0].header
                self.set_base_headers(hdr)
                # rectification and wavelength calibration
                reduced_mos_image = apply_rectwv_coeff(
                    reduced_image,
                    list_rectwv_coeff[i]
                )
                self.save_intermediate_img(
                    reduced_mos_image,
                    'reduced_mos_image_' + char + '_' + frame.filename[:10] +
                    '.fits')

                # get useful pixels in the wavelength direction
                if char == 'A' and first_a:
                    xisok_a = useful_mos_xpixels(
                        reduced_mos_image=reduced_mos_image,
                        vertical_pix_region=vertical_pix_region_a,
                        npix_removed_around_ohlines=npix_removed_around_ohlines,
                        list_valid_wvregions=list_valid_wvregions_a,
                        debugplot=0
                    )
                elif char == 'B' and first_b:
                    xisok_b = useful_mos_xpixels(
                        reduced_mos_image=reduced_mos_image,
                        vertical_pix_region=vertical_pix_region_b,
                        npix_removed_around_ohlines=npix_removed_around_ohlines,
                        list_valid_wvregions=list_valid_wvregions_b,
                        debugplot=0
                    )

                if char == 'A':
                    nsmin = vertical_pix_region_a[0]
                    nsmax = vertical_pix_region_a[1]
                    xisok = xisok_a
                elif char == 'B':
                    nsmin = vertical_pix_region_b[0]
                    nsmax = vertical_pix_region_b[1]
                    xisok = xisok_b
                else:
                    raise ValueError('Unexpected char value: {}'.format(char))

                # initial slitlet region
                slitlet2d = reduced_mos_image[0].data[(nsmin - 1):nsmax, :].copy()

                # selected wavelength regions after blocking OH lines
                slitlet2d_blocked = slitlet2d[:, xisok]

                # apply median filter in the X direction
                if nwidth_medfilt > 1:
                    slitlet2d_blocked_smoothed = ndimage.filters.median_filter(
                        slitlet2d_blocked, size=(1, nwidth_medfilt)
                    )
                else:
                    slitlet2d_blocked_smoothed = slitlet2d_blocked

                # ToDo: choose best option
                # convert to 1D series of spatial profiles
                # spatial_prof = slitlet2d_blocked_smoothed.transpose().ravel()
                profile = np.sum(slitlet2d_blocked_smoothed, axis=1)
                if char == 'A' and first_a:
                    reference_profile_a = profile.copy()
                elif char== 'B' and first_b:
                    reference_profile_b = profile.copy()
                # ToDo: save profile using pickle
                # para leer
                # with open('test.pkl', 'rb') as f:
                #     x = pickle.load(f)
                import pickle
                with open(frame.filename[:10] + '_profile.pkl', 'wb') as f:
                    pickle.dump(profile, f)

                # crosscorrelation to find offset
                naround_zero = (nsmax - nsmin) // 2
                if char == 'A':
                    reference_profile = reference_profile_a
                elif char == 'B':
                    reference_profile = reference_profile_b
                else:
                    raise ValueError('Unexpected char value: {}'.format(char))
                offset, fpeak = periodic_corr1d(
                    reference_profile,
                    profile,
                    fminmax=None,
                    zero_padding=0,
                    nfit_peak=7,
                    naround_zero=naround_zero,
                    plottitle='Image {} (type {})'.format(i + 1, char),
                    debugplot=22
                )
                # round to 4 decimal places
                if abs(offset) < 1E-4:
                    offset = 0.0   # avoid -0.0
                else:
                    offset = round(offset, 4)
                list_offsets.append(offset)

                # end of loop
                if char == 'A' and first_a:
                    first_a = False
                elif char == 'B' and first_b:
                    first_b = False
            else:
                list_offsets.append(0.0)

        list_a = []
        list_b = []
        for i, (char, offset) in enumerate(zip(full_set, list_offsets)):
            frame = rinput.obresult.frames[i]
            self.logger.info('image {} ({} of {})'.format(char, i + 1,
                                                          nimages))
            self.logger.info('image: {}'.format(frame.filename))
            with frame.open() as f:
                base_header = f[0].header.copy()
                data = f[0].data.astype('double')
            if int(offset*1000 + 0.5) != 0:
                self.logger.info(
                    'correcting vertical offset (pixesl): {}'.format(offset))
                data = shift_image2d(data, yoffset=-offset).astype('float32')
            hdu = fits.PrimaryHDU(data, header=base_header)
            hdu.header['UUID'] = str(uuid.uuid1())
            hdul = fits.HDUList([hdu])
            for hdu in f[1:]:
                hdul.append(hdu.copy())
            # basic reduction
            reduced_image = flow(hdul)
            hdr = reduced_image[0].header
            self.set_base_headers(hdr)
            # rectification and wavelength calibration
            reduced_mos_image = apply_rectwv_coeff(
                reduced_image,
                list_rectwv_coeff[i]
            )
            self.save_intermediate_img(
                reduced_mos_image,
                'reduced_mos_image_refined_' + char + '_' +
                frame.filename[:10] + '.fits')
            # store reduced_mos_image
            if char == 'A':
                list_a.append(reduced_mos_image)
            elif char == 'B':
                list_b.append(reduced_mos_image)
            else:
                raise ValueError('Unexpected char value: {}'.format(char))

        # combination method
        fmethod = getattr(combine, rinput.method)

        # final combination of A images
        reduced_mos_image_a = combine_imgs(
            list_a,
            method=fmethod,
            method_kwargs=rinput.method_kwargs,
            errors=False,
            prolog=None
        )
        self.save_intermediate_img(reduced_mos_image_a,
                                   'reduced_mos_image_a.fits')

        # final combination of B images
        reduced_mos_image_b = combine_imgs(
            list_b,
            method=fmethod,
            method_kwargs=rinput.method_kwargs,
            errors=False,
            prolog=None
        )
        self.save_intermediate_img(reduced_mos_image_b,
                                   'reduced_mos_image_b.fits')

        header_a = reduced_mos_image_a[0].header
        header_b = reduced_mos_image_b[0].header
        data_a = reduced_mos_image_a[0].data.astype('float32')
        data_b = reduced_mos_image_b[0].data.astype('float32')
        reduced_mos_abba_data = data_a - data_b

        # update reduced image header
        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(fname.open()) for fname in
                     rinput.obresult.frames]
            reduced_mos_abba = self.create_mos_abba_image(
                hduls,
                reduced_mos_abba_data,
                header_a, header_b,
                pattern,
                full_set,
                list_offsets
            )

        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(list_rectwv_coeff[0])
            save_spectral_lines_ds9(list_rectwv_coeff[0])

        # ToDo: estas imagenes no tienen la calibracion en longitud de onda
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

    def create_mos_abba_image(self, hduls, reduced_mos_abba_data,
                              header_a, header_b,
                              pattern, full_set, list_offsets):
        # Copy header of first image
        base_header = hduls[0][0].header.copy()

        hdu = fits.PrimaryHDU(reduced_mos_abba_data, header=base_header)
        self.set_base_headers(hdu.header)

        self.logger.debug('update result header')
        crpix1 = header_a['crpix1']
        crval1 = header_a['crval1']
        cdelt1 = header_a['cdelt1']
        if crpix1 != header_b['crpix1'] or \
            crval1 != header_b['crval1'] or cdelt1 != header_b['cdelt1']:
            raise ValueError('Unexpected differences in wavelength '
                             'calibration  parameters')

        # update wavelength calibration in FITS header
        for keyword in ['crval1', 'crpix1', 'crval2', 'crpix2']:
            if keyword in hdu.header:
                hdu.header.remove(keyword)
        hdu.header['crpix1'] = (crpix1, 'reference pixel')
        hdu.header['crval1'] = (crval1, 'central wavelength at crpix1')
        hdu.header['cdelt1'] = (cdelt1, 'linear dispersion (Angstrom/pixel)')
        hdu.header['cunit1'] = ('Angstrom', 'units along axis1')
        hdu.header['ctype1'] = 'WAVELENGTH'
        hdu.header['crpix2'] = (0.0, 'reference pixel')
        hdu.header['crval2'] = (0.0, 'central value at crpix2')
        hdu.header['cdelt2'] = (1.0, 'increment')
        hdu.header['ctype2'] = 'PIXEL'
        hdu.header['cunit2'] = ('Pixel', 'units along axis2')
        for keyword in ['cd1_1', 'cd1_2', 'cd2_1', 'cd2_2',
                        'PCD1_1', 'PCD1_2', 'PCD2_1', 'PCD2_2',
                        'PCRPIX1', 'PCRPIX2']:
            if keyword in hdu.header:
                hdu.header.remove(keyword)

        # update additional keywords
        hdu.header['UUID'] = str(uuid.uuid1())
        hdu.header['OBSMODE'] = pattern + ' pattern'
        hdu.header['TSUTC2'] = hduls[-1][0].header['TSUTC2']
        hdu.header['NUM-NCOM'] = (len(hduls), 'Number of combined frames')

        # update history
        hdu.header['HISTORY'] = "Processed " + pattern + " pattern"
        hdu.header['HISTORY'] = '--- Reduction of A images ---'
        for line in header_a['HISTORY']:
            hdu.header['HISTORY'] = line
        hdu.header['HISTORY'] = '--- Reduction of B images ---'
        for line in header_b['HISTORY']:
            hdu.header['HISTORY'] = line
        hdu.header['HISTORY'] = '--- Combination of ABBA images ---'
        dm = emirdrp.datamodel.EmirDataModel()
        for img, key, offset in zip(hduls, full_set, list_offsets):
            imgid = dm.get_imgid(img)
            hdu.header['HISTORY'] = \
                "Image '{}' is '{}', with voffset {}".format(imgid, key, offset)
        result = fits.HDUList([hdu])
        return result

    def set_base_headers(self, hdr):
        newhdr = super(ABBASpectraRectwv, self).set_base_headers(hdr)
        # Update EXP to 0
        newhdr['EXP'] = 0
        return newhdr

    def extract_tags_from_obsres(self, obsres, tag_keys):
        # this function is necessary to make use of recipe requirements
        # that are provided as lists
        final_tags = []
        for frame in obsres.frames:
            ref_img = frame.open()
            tag = self.extract_tags_from_ref(
                ref_img,
                tag_keys,
                base=obsres.labels
            )
            final_tags.append(tag)
        return final_tags


class ABBASpectraFastRectwv(EmirRecipe):
    """Process AB or ABBA images applying a single wavelength calibration

    Note that the wavelength calibration has already been determined.

    """

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    rectwv_coeff = reqs.RectWaveCoeffRequirement()

    pattern = Parameter(
        'ABBA',
        description='Observation pattern',
        choices=['AB', 'ABBA']
    )
    # do not allow 'sum' as combination method in order to avoid
    # confusion with number of counts in resulting image
    method = Parameter(
        'mean',
        description='Combination method',
        choices=['mean', 'median', 'sigmaclip']
    )
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method'
    )

    reduced_mos_abba = Result(prods.ProcessedMOS)

    def run(self, rinput):
        nimages = len(rinput.obresult.frames)
        pattern = rinput.pattern
        pattern_length = len(pattern)
        # check pattern sequence matches number of images
        if nimages % pattern_length != 0:
            raise ValueError('Number of images is not a multiple of pattern '
                             'length: {}, {}'.format(nimages, pattern_length))
        nsequences = nimages // pattern_length

        rectwv_coeff = rinput.rectwv_coeff
        self.logger.info(rectwv_coeff)
        self.logger.info('observation pattern: {}'.format(pattern))
        self.logger.info('nsequences.........: {}'.format(nsequences))

        full_set = pattern * nsequences
        self.logger.info('full set of images.: {}'.format(full_set))

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
                method_kwargs=rinput.method_kwargs,
                errors=False,
                prolog=None
            )
        reduced_image_a = flow(reduced_image_a)
        hdr = reduced_image_a[0].header
        self.set_base_headers(hdr)

        # basic reduction of B images
        list_b = [rinput.obresult.frames[i] for i, char in
                  enumerate(full_set) if char == 'B']
        with contextlib.ExitStack() as stack:
            self.logger.info('starting basic reduction of B images')
            hduls = [stack.enter_context(fname.open()) for fname in list_b]
            reduced_image_b = combine_imgs(
                hduls,
                method=fmethod,
                method_kwargs=rinput.method_kwargs,
                errors=False,
                prolog=None
            )
        reduced_image_b = flow(reduced_image_b)
        hdr = reduced_image_b[0].header
        self.set_base_headers(hdr)

        # save intermediate reduced_image_a and reduced_image_b
        self.save_intermediate_img(reduced_image_a, 'reduced_image_a.fits')
        self.save_intermediate_img(reduced_image_b, 'reduced_image_b.fits')

        # computation of A-B
        data_a = reduced_image_a[0].data.astype('float32')
        data_b = reduced_image_b[0].data.astype('float32')
        reduced_data = data_a - data_b

        # update reduced image header
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
        self.logger.info('begin rect.+wavecal. reduction of ABBA spectra')
        reduced_mos_abba = apply_rectwv_coeff(
            reduced_image,
            rectwv_coeff
        )

        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(rectwv_coeff)
            save_spectral_lines_ds9(rectwv_coeff)

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
        newhdr = super(ABBASpectraFastRectwv, self).set_base_headers(hdr)
        # Update EXP to 0
        newhdr['EXP'] = 0
        return newhdr

    def extract_tags_from_obsres(self, obsres, tag_keys):
        # this function is necessary to make use of recipe requirements
        # that are provided as lists
        final_tags = []
        for frame in obsres.frames:
            ref_img = frame.open()
            tag = self.extract_tags_from_ref(
                ref_img,
                tag_keys,
                base=obsres.labels
            )
            final_tags.append(tag)
        return final_tags
