#
# Copyright 2019-2020 Universidad Complutense de Madrid
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
import logging
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
    """Compute offsets and spatial scales from WCS info in image headers.

    Compute offset (arcsec) between the first image in the sequence and
    the rest of the images, as well as the spatial scale (arcsec/pix)
    within each individual image.

    """

    nimages = len(hduls)
    c0 = None
    list_sep_arcsec = []
    list_spatial_scales_arcsecperpix = []
    for i in range(nimages):
        wcsh = wcs.WCS(hduls[i][0].header)
        ra, dec = wcsh.wcs_pix2world(EMIR_NAXIS1 / 2 + 0.5,
                                     EMIR_NAXIS2 / 2 + 0.5, 1)
        ra_, dec_ = wcsh.wcs_pix2world(EMIR_NAXIS1 / 2 + 0.5,
                                       EMIR_NAXIS2 / 2 + 1.5, 1)
        if i == 0:
            c0 = SkyCoord(ra, dec, unit="deg")
        ci = SkyCoord(ra, dec, unit="deg")
        ci_ = SkyCoord(ra_, dec_, unit="deg")
        sep = ci.separation(c0).arcsec
        scale = ci.separation(ci_).arcsec
        list_sep_arcsec.append(round(sep, 4))
        list_spatial_scales_arcsecperpix.append(round(scale, 7))

    sep_arcsec = np.array(list_sep_arcsec)
    spatial_scales_arcsecperpix = np.array(list_spatial_scales_arcsecperpix)

    return sep_arcsec, spatial_scales_arcsecperpix


class ABBASpectraRectwv(EmirRecipe):
    """Process images in AB or ABBA  mode applying wavelength calibration

    Note that in this case the wavelength calibration has already been
    determined.

    """

    logger = logging.getLogger(__name__)

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
    reduced_mos_abba_combined = Result(prods.ProcessedMOS)

    def run(self, rinput):

        nimages = len(rinput.obresult.frames)
        pattern = rinput.pattern
        pattern_length = len(pattern)

        # check combination method
        if rinput.method != 'sigmaclip':
            if rinput.method_kwargs != {}:
                raise ValueError('Unexpected method_kwargs={}'.format(
                    rinput.method_kwargs))

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
                        raise ValueError(
                            'list_rectwv_coeff contains coefficients for '
                            'different {}s'.format(item)
                        )
        else:
            raise ValueError("Unexpected error!")

        # grism and filter names
        grism_name = list_rectwv_coeff[0].tags['grism']
        filter_name = list_rectwv_coeff[0].tags['filter']

        # compute offsets from WCS info in image headers
        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(fname.open()) for fname in
                     rinput.obresult.frames]
            sep_arcsec, spatial_scales = compute_wcs_offsets(hduls)

        sep_pixel = np.round(sep_arcsec / spatial_scales, 6)

        for i in range(nimages):
            self.logger.info(list_rectwv_coeff[i])
        self.logger.info('observation pattern..................: {}'.format(
            pattern))
        self.logger.info('nsequences...........................: {}'.format(
            nsequences))
        self.logger.info('offsets (arcsec, from WCS)...........: {}'.format(
            sep_arcsec))
        self.logger.info('spatial scales (arcsec/pix, from WCS): {}'.format(
            spatial_scales))
        self.logger.info('spatial scales (pixels, from WCS)....: {}'.format(
            sep_pixel))

        full_set = pattern * nsequences
        self.logger.info('full set of images...................: {}'.format(
            full_set))

        # basic parameters to determine useful pixels in the wavelength
        # direction
        dict_rtas = rinput.refine_target_along_slitlet

        valid_keys = [
            'npix_removed_near_ohlines',
            'nwidth_medfilt',
            'save_individual_images',
            'ab_different_target',
            'vpix_region_a_target',
            'vpix_region_a_sky',
            'vpix_region_b_target',
            'vpix_region_b_sky',
            'list_valid_wvregions_a',
            'list_valid_wvregions_b'
        ]
        for dumkey in dict_rtas.keys():
            if dumkey not in valid_keys:
                raise ValueError('Unexpected key={}'.format(dumkey))

        if 'vpix_region_a_target' in dict_rtas.keys():
            vpix_region_a_target = dict_rtas['vpix_region_a_target']
        else:
            vpix_region_a_target = None

        if 'vpix_region_a_sky' in dict_rtas.keys():
            vpix_region_a_sky = dict_rtas['vpix_region_a_sky']
        else:
            vpix_region_a_sky = None

        if 'vpix_region_b_target' in dict_rtas.keys():
            vpix_region_b_target = dict_rtas['vpix_region_b_target']
        else:
            vpix_region_b_target = None

        if 'vpix_region_b_sky' in dict_rtas.keys():
            vpix_region_b_sky = dict_rtas['vpix_region_b_sky']
        else:
            vpix_region_b_sky = None

        if 'ab_different_target' in dict_rtas.keys():
            ab_different_target = int(dict_rtas['ab_different_target'])
            if ab_different_target not in [-1, 0, 1, 9]:
                raise ValueError('Invalid ab_different_target={} value'.format(
                    ab_different_target))
        else:
            raise ValueError('Missing ab_different_target value')

        try:
            npix_removed_near_ohlines = \
                int(dict_rtas['npix_removed_near_ohlines'])
        except KeyError:
            npix_removed_near_ohlines = 0
        except ValueError:
            raise ValueError(
                'wrong value: npix_removed_near_ohlines='
                '{}'.format(dict_rtas['npix_removed_near_ohlines'])
            )

        if 'list_valid_wvregions_a' in dict_rtas.keys():
            if vpix_region_a_target is None:
                raise ValueError('Unexpected list_valid_wvregions_a when '
                                 'vpix_region_a_target is not set')
            list_valid_wvregions_a = dict_rtas['list_valid_wvregions_a']
        else:
            list_valid_wvregions_a = None

        if 'list_valid_wvregions_b' in dict_rtas.keys():
            if vpix_region_b_target is None:
                raise ValueError('Unexpected list_valid_wvregions_b when '
                                 'vpix_region_b_target is not set')

            list_valid_wvregions_b = dict_rtas['list_valid_wvregions_b']
        else:
            list_valid_wvregions_b = None

        try:
            nwidth_medfilt = int(dict_rtas['nwidth_medfilt'])
        except KeyError:
            nwidth_medfilt = 0
        except ValueError:
            raise ValueError('wrong value: nwidth_medfilt={}'.format(
                dict_rtas['nwidth_medfilt']))

        try:
            save_individual_images = int(dict_rtas['save_individual_images'])
        except KeyError:
            save_individual_images = 0
        except ValueError:
            raise ValueError('wrong value: save_individual_images={}'.format(
                dict_rtas['save_individual_images']))

        self.logger.info('npix_removed_near_ohlines: {}'.format(
            npix_removed_near_ohlines))
        self.logger.info('nwidth_medfilt: {}'.format(nwidth_medfilt))
        self.logger.info('vpix_region_a_target: {}'.format(
            vpix_region_a_target))
        self.logger.info('vpix_region_b_target: {}'.format(
            vpix_region_b_target))
        self.logger.info('list_valid_wvregions_a: {}'.format(
            list_valid_wvregions_a))
        self.logger.info('list_valid_wvregions_b: {}'.format(
            list_valid_wvregions_b))

        # build object to proceed with bpm, bias, dark and flat
        flow = self.init_filters(rinput)

        # basic reduction, rectification and wavelength calibration of
        # all the individual images
        list_reduced_mos_images = []
        self.logger.info('starting reduction of individual images')
        for i, char in enumerate(full_set):
            frame = rinput.obresult.frames[i]
            self.logger.info('image {} ({} of {})'.format(char, i + 1,
                                                          nimages))
            self.logger.info('image: {}'.format(frame.filename))
            with frame.open() as f:
                newimg = fits.HDUList([ext.copy() for ext in f])
            base_header = newimg[0].header
            grism_name_ = base_header['grism']
            if grism_name_ != grism_name:
                raise ValueError('Incompatible grism name in '
                                 'rectwv_coeff.json file and FITS image')
            filter_name_ = base_header['filter']
            if filter_name_ != filter_name:
                raise ValueError('Incompatible filter name in '
                                 'rectwv_coeff.json file and FITS image')
            hdu = newimg[0]
            hdu.header['UUID'] = str(uuid.uuid1())
            # basic reduction
            reduced_image = flow(newimg)
            hdr = reduced_image[0].header
            self.set_base_headers(hdr)
            # rectification and wavelength calibration
            reduced_mos_image = apply_rectwv_coeff(
                reduced_image,
                list_rectwv_coeff[i]
            )
            if save_individual_images != 0:
                self.save_intermediate_img(
                    reduced_mos_image,
                    'reduced_mos_image_' + char + '_' +
                    frame.filename[:10] + '.fits'
                )
            list_reduced_mos_images.append(reduced_mos_image)

        # intermediate PDF file with crosscorrelation plots
        if self.intermediate_results:
            from matplotlib.backends.backend_pdf import PdfPages
            pdf = PdfPages('crosscorrelation_ab.pdf')
        else:
            pdf = None

        # compute offsets between images
        self.logger.info('computing offsets between individual images')
        first_a = True
        first_b = True
        xisok_a = None
        xisok_b = None
        reference_profile_a = None
        reference_profile_b = None
        refine_a = vpix_region_a_target is not None
        refine_b = vpix_region_b_target is not None
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
                data = list_reduced_mos_images[i][0].data.copy()
                data_sky = list_reduced_mos_images[isky][0].data
                data -= data_sky
                base_header = list_reduced_mos_images[i][0].header

                # get useful pixels in the wavelength direction
                if char == 'A' and first_a:
                    xisok_a = useful_mos_xpixels(
                        data,
                        base_header,
                        vpix_region=vpix_region_a_target,
                        npix_removed_near_ohlines=npix_removed_near_ohlines,
                        list_valid_wvregions=list_valid_wvregions_a,
                        debugplot=0
                    )
                elif char == 'B' and first_b:
                    xisok_b = useful_mos_xpixels(
                        data,
                        base_header,
                        vpix_region=vpix_region_b_target,
                        npix_removed_near_ohlines=npix_removed_near_ohlines,
                        list_valid_wvregions=list_valid_wvregions_b,
                        debugplot=0
                    )

                if char == 'A':
                    nsmin = vpix_region_a_target[0]
                    nsmax = vpix_region_a_target[1]
                    xisok = xisok_a
                    if vpix_region_a_sky is not None:
                        nsmin_sky = vpix_region_a_sky[0]
                        nsmax_sky = vpix_region_a_sky[1]
                        skysubtraction = True
                    else:
                        nsmin_sky = None
                        nsmax_sky = None
                        skysubtraction = False
                elif char == 'B':
                    nsmin = vpix_region_b_target[0]
                    nsmax = vpix_region_b_target[1]
                    xisok = xisok_b
                    if vpix_region_b_sky is not None:
                        nsmin_sky = vpix_region_b_sky[0]
                        nsmax_sky = vpix_region_b_sky[1]
                        skysubtraction = True
                    else:
                        nsmin_sky = None
                        nsmax_sky = None
                        skysubtraction = False
                else:
                    raise ValueError('Unexpected char value: {}'.format(char))

                # initial slitlet region
                slitlet2d = data[(nsmin - 1):nsmax, :].copy()
                if skysubtraction:
                    slitlet2d_sky = data[(nsmin_sky - 1):nsmax_sky, :].copy()
                    median_sky = np.median(slitlet2d_sky, axis=0)
                    slitlet2d -= median_sky

                # selected wavelength regions after blocking OH lines
                slitlet2d_blocked = slitlet2d[:, xisok]

                # apply median filter in the X direction
                if nwidth_medfilt > 1:
                    slitlet2d_blocked_smoothed = ndimage.filters.median_filter(
                        slitlet2d_blocked, size=(1, nwidth_medfilt)
                    )
                else:
                    slitlet2d_blocked_smoothed = slitlet2d_blocked

                # ---
                # convert to 1D series of spatial profiles (note: this
                # option does not work very well when the signal is low;
                # a simple mean profile does a better job)
                # profile = slitlet2d_blocked_smoothed.transpose().ravel()
                # ---
                profile = np.mean(slitlet2d_blocked_smoothed, axis=1)
                if char == 'A' and first_a:
                    reference_profile_a = profile.copy()
                elif char == 'B' and first_b:
                    reference_profile_b = profile.copy()
                # ---
                # # To write pickle file
                # import pickle
                # with open(frame.filename[:10] + '_profile.pkl', 'wb') as f:
                #     pickle.dump(profile, f)
                # # To read pickle file
                # with open('test.pkl', 'rb') as f:
                #     x = pickle.load(f)
                # ---

                # crosscorrelation to find offset
                naround_zero = (nsmax - nsmin) // 3
                if char == 'A':
                    reference_profile = reference_profile_a
                elif char == 'B':
                    reference_profile = reference_profile_b
                else:
                    raise ValueError('Unexpected char value: {}'.format(char))
                offset, fpeak = periodic_corr1d(
                    sp_reference=reference_profile,
                    sp_offset=profile,
                    remove_mean=False,
                    frac_cosbell=0.10,
                    zero_padding=11,
                    fminmax=None,
                    nfit_peak=5,
                    naround_zero=naround_zero,
                    sp_label='spatial profile',
                    plottitle='Image #{} (type {}), {}'.format(
                        i + 1, char, frame.filename[:10]),
                    pdf=pdf
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
        self.logger.info('computed offsets: {}'.format(list_offsets))

        self.logger.info('correcting vertical offsets between individual '
                         'images')
        list_a = []
        list_b = []
        for i, (char, offset) in enumerate(zip(full_set, list_offsets)):
            frame = rinput.obresult.frames[i]
            self.logger.info('image {} ({} of {})'.format(char, i + 1,
                                                          nimages))
            self.logger.info('image: {}'.format(frame.filename))
            reduced_mos_image = list_reduced_mos_images[i]
            data = reduced_mos_image[0].data
            base_header = reduced_mos_image[0].header
            self.logger.info(
                'correcting vertical offset (pixesl): {}'.format(offset))
            if offset != 0:
                reduced_mos_image[0].data = shift_image2d(
                    data,
                    yoffset=-offset
                ).astype('float32')
            base_header['HISTORY'] = 'Applying voffset_pix {}'.format(offset)
            if save_individual_images != 0:
                self.save_intermediate_img(
                    reduced_mos_image,
                    'reduced_mos_image_refined_' + char + '_' +
                    frame.filename[:10] + '.fits'
                )

            # store reduced_mos_image
            if char == 'A':
                list_a.append(reduced_mos_image)
            elif char == 'B':
                list_b.append(reduced_mos_image)
            else:
                raise ValueError('Unexpected char value: {}'.format(char))

        # combination method
        method = getattr(combine, rinput.method)
        method_kwargs = rinput.method_kwargs

        # final combination of A images
        self.logger.info('combining individual A images')
        reduced_mos_image_a = combine_imgs(
            list_a,
            method=method,
            method_kwargs=method_kwargs,
            errors=False,
            prolog=None
        )
        self.save_intermediate_img(reduced_mos_image_a,
                                   'reduced_mos_image_a.fits')

        # final combination of B images
        self.logger.info('combining individual B images')
        reduced_mos_image_b = combine_imgs(
            list_b,
            method=method,
            method_kwargs=method_kwargs,
            errors=False,
            prolog=None
        )
        self.save_intermediate_img(reduced_mos_image_b,
                                   'reduced_mos_image_b.fits')

        self.logger.info('mixing A and B spectra')
        header_a = reduced_mos_image_a[0].header
        header_b = reduced_mos_image_b[0].header
        data_a = reduced_mos_image_a[0].data.astype('float32')
        data_b = reduced_mos_image_b[0].data.astype('float32')

        reduced_mos_abba_data = data_a - data_b

        # update reduced mos image header
        reduced_mos_abba = self.create_mos_abba_image(
            rinput,
            dict_rtas,
            reduced_mos_abba_data,
            header_a, header_b,
            pattern,
            full_set,
            list_offsets,
            voffset_pix=0
        )

        # combine A and B data by shifting B on top of A
        if abs(ab_different_target) == 0:
            len_prof_a = len(reference_profile_a)
            len_prof_b = len(reference_profile_b)
            if len_prof_a == len_prof_b:
                reference_profile = reference_profile_a
                profile = reference_profile_b
                naround_zero = len_prof_a // 3
            elif len_prof_a > len_prof_b:
                ndiff = len_prof_a - len_prof_b
                reference_profile = reference_profile_a
                profile = np.concatenate(
                    (reference_profile_b, np.zeros(ndiff, dtype='float'))
                )
                naround_zero = len_prof_a // 3
            else:
                ndiff = len_prof_b - len_prof_a
                reference_profile = np.concatenate(
                    (reference_profile_a, np.zeros(ndiff, dtype='float'))
                )
                profile = reference_profile_b
                naround_zero = len_prof_b // 3
            offset, fpeak = periodic_corr1d(
                sp_reference=reference_profile,
                sp_offset=profile,
                remove_mean=False,
                frac_cosbell=0.10,
                zero_padding=11,
                fminmax=None,
                nfit_peak=5,
                naround_zero=naround_zero,
                sp_label='spatial profile',
                plottitle='Comparison of A and B profiles',
                pdf=pdf
            )
            voffset_pix = vpix_region_b_target[0] - vpix_region_a_target[0]
            voffset_pix += offset
        elif abs(ab_different_target) == 1:
            # apply nominal offset (from WCS info) between first A
            # and first B
            voffset_pix = sep_pixel[1] * ab_different_target
        elif ab_different_target == 9:
            voffset_pix = None
        else:
            raise ValueError('Invalid ab_different_target={}'.format(
                ab_different_target))

        # close output PDF file
        if pdf is not None:
            pdf.close()

        if voffset_pix is not None:
            self.logger.info(
                'correcting vertical offset (pixesl): {}'.format(voffset_pix))
            shifted_a_minus_b_data = shift_image2d(
                reduced_mos_abba_data,
                yoffset=-voffset_pix,
            ).astype('float32')
            reduced_mos_abba_combined_data = \
                reduced_mos_abba_data - shifted_a_minus_b_data
            # scale signal to exposure of a single image
            reduced_mos_abba_combined_data /= 2.0
        else:
            reduced_mos_abba_combined_data = None

        # update reduced mos combined image header
        reduced_mos_abba_combined = self.create_mos_abba_image(
            rinput,
            dict_rtas,
            reduced_mos_abba_combined_data,
            header_a, header_b,
            pattern,
            full_set,
            list_offsets,
            voffset_pix
        )

        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(list_rectwv_coeff[0])
            save_spectral_lines_ds9(list_rectwv_coeff[0])

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
        result = self.create_result(
            reduced_mos_abba=reduced_mos_abba,
            reduced_mos_abba_combined=reduced_mos_abba_combined
        )
        return result

    def create_mos_abba_image(self, rinput, dict_rtas, reduced_mos_abba_data,
                              header_a, header_b,
                              pattern, full_set, list_offsets, voffset_pix):
        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(fname.open()) for fname in
                     rinput.obresult.frames]
            # Copy the first image
            result_img = fits.HDUList([ext.copy() for ext in hduls[0]])
            hdu = result_img[0]
            hdu.data = reduced_mos_abba_data
            self.set_base_headers(hdu.header)

            # check consistency of wavelength calibration paramenters
            for param in ['crpix1', 'crval1', 'cdelt1']:
                if header_a[param] != header_b[param]:
                    raise ValueError('Headers of A and B images have different '
                                     'values of {}'.format(param))
            self.logger.debug('update result header')
            crpix1 = header_a['crpix1']
            crval1 = header_a['crval1']
            cdelt1 = header_a['cdelt1']

            # update wavelength calibration in FITS header
            for keyword in ['crval1', 'crpix1', 'crval2', 'crpix2']:
                if keyword in hdu.header:
                    hdu.header.remove(keyword)
            hdu.header['crpix1'] = (crpix1, 'reference pixel')
            hdu.header['crval1'] = (crval1, 'central wavelength at crpix1')
            hdu.header['cdelt1'] = (cdelt1, 'linear dispersion (Angstrom/pixel)')
            hdu.header['cunit1'] = ('Angstrom', 'units along axis1')
            hdu.header['ctype1'] = 'WAVE'
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
            for key in dict_rtas:
                hdu.header['HISTORY'] = '{}: {}'.format(key, dict_rtas[key])
            dm = emirdrp.datamodel.EmirDataModel()
            for img, key, offset in zip(hduls, full_set, list_offsets):
                imgid = dm.get_imgid(img)
                hdu.header['HISTORY'] = \
                    "Image '{}' is '{}', with voffset_pix {}".format(
                        imgid, key, offset)

        if voffset_pix is not None and voffset_pix != 0:
            hdu.header['HISTORY'] = '--- Combination of AB spectra ---'
            hdu.header['HISTORY'] = "voffset_pix between A and B {}".format(
                voffset_pix)
        hdu.header.add_history('--- rinput.stored() (BEGIN) ---')
        for item in rinput.stored():
            value = getattr(rinput, item)
            cline = '{}: {}'.format(item, value)
            hdu.header.add_history(cline)
        hdu.header.add_history('--- rinput.stored() (END) ---')
        return result_img

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

    logger = logging.getLogger(__name__)

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
    voffset_pix = Parameter(
        0.0,
        description='Shift (pixels) to move A into B',
        optional=True
    )

    reduced_mos_abba = Result(prods.ProcessedMOS)
    reduced_mos_abba_combined = Result(prods.ProcessedMOS)

    def run(self, rinput):

        nimages = len(rinput.obresult.frames)
        pattern = rinput.pattern
        pattern_length = len(pattern)

        # check combination method
        if rinput.method != 'sigmaclip':
            if rinput.method_kwargs != {}:
                raise ValueError('Unexpected method_kwargs={}'.format(
                    rinput.method_kwargs))

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
        method = getattr(combine, rinput.method)
        method_kwargs = rinput.method_kwargs

        # basic reduction of A images
        list_a = [rinput.obresult.frames[i] for i, char in enumerate(full_set)
                  if char == 'A']
        with contextlib.ExitStack() as stack:
            self.logger.info('starting basic reduction of A images')
            hduls = [stack.enter_context(fname.open()) for fname in list_a]
            reduced_image_a = combine_imgs(
                hduls,
                method=method,
                method_kwargs=method_kwargs,
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
                method=method,
                method_kwargs=method_kwargs,
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
        header_a = reduced_image_a[0].header
        header_b = reduced_image_b[0].header
        data_a = reduced_image_a[0].data.astype('float32')
        data_b = reduced_image_b[0].data.astype('float32')
        reduced_data = data_a - data_b

        # update reduced image header
        reduced_image = self.create_reduced_image(
            rinput,
            reduced_data,
            header_a,
            header_b,
            rinput.pattern,
            full_set,
            voffset_pix=0,
            header_mos_abba=None,
        )

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # apply rectification and wavelength calibration
        self.logger.info('begin rect.+wavecal. reduction of ABBA spectra')
        reduced_mos_abba = apply_rectwv_coeff(
            reduced_image,
            rectwv_coeff
        )
        header_mos_abba = reduced_mos_abba[0].header

        # combine A and B data by shifting B on top of A
        voffset_pix = rinput.voffset_pix

        if voffset_pix is not None and voffset_pix != 0:
            self.logger.info(
                'correcting vertical offset (pixesl): {}'.format(voffset_pix))
            reduced_mos_abba_data = reduced_mos_abba[0].data.astype('float32')
            shifted_a_minus_b_data = shift_image2d(
                reduced_mos_abba_data,
                yoffset=-voffset_pix,
            ).astype('float32')
            reduced_mos_abba_combined_data = \
                reduced_mos_abba_data - shifted_a_minus_b_data
            # scale signal to exposure of a single image
            reduced_mos_abba_combined_data /= 2.0
        else:
            reduced_mos_abba_combined_data = None

        # update reduced combined image header
        reduced_mos_abba_combined = self.create_reduced_image(
            rinput,
            reduced_mos_abba_combined_data,
            header_a,
            header_b,
            rinput.pattern,
            full_set,
            voffset_pix,
            header_mos_abba=header_mos_abba
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
        result = self.create_result(
            reduced_mos_abba=reduced_mos_abba,
            reduced_mos_abba_combined=reduced_mos_abba_combined
        )
        return result

    def create_reduced_image(self, rinput, reduced_data,
                             header_a, header_b,
                             pattern, full_set,
                             voffset_pix, header_mos_abba):
        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(fname.open()) for fname in
                     rinput.obresult.frames]

            # Copy the first image
            result_img = fits.HDUList([ext.copy() for ext in hduls[0]])
            hdu = result_img[0]
            hdu.data = reduced_data
            self.set_base_headers(hdu.header)

            self.logger.debug('update result header')
            if header_mos_abba is not None:
                self.logger.debug('update result header')
                crpix1 = header_mos_abba['crpix1']
                crval1 = header_mos_abba['crval1']
                cdelt1 = header_mos_abba['cdelt1']

                # update wavelength calibration in FITS header
                for keyword in ['crval1', 'crpix1', 'crval2', 'crpix2']:
                    if keyword in hdu.header:
                        hdu.header.remove(keyword)
                hdu.header['crpix1'] = (crpix1, 'reference pixel')
                hdu.header['crval1'] = (crval1, 'central wavelength at crpix1')
                hdu.header['cdelt1'] = \
                    (cdelt1, 'linear dispersion (Angstrom/pixel)')
                hdu.header['cunit1'] = ('Angstrom', 'units along axis1')
                hdu.header['ctype1'] = 'WAVE'
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
            hdu.header['history'] = "Processed " + pattern + " pattern"
            hdu.header['NUM-NCOM'] = (len(hduls), 'Number of combined frames')

            # update history
            dm = emirdrp.datamodel.EmirDataModel()
            for img, key in zip(hduls, full_set):
                imgid = dm.get_imgid(img)
                hdu.header['HISTORY'] = "Image '{}' is '{}'".format(imgid, key)

        hdu.header['HISTORY'] = "Processed " + pattern + " pattern"
        hdu.header['HISTORY'] = '--- Reduction of A images ---'
        for line in header_a['HISTORY']:
            hdu.header['HISTORY'] = line
        hdu.header['HISTORY'] = '--- Reduction of B images ---'
        for line in header_b['HISTORY']:
            hdu.header['HISTORY'] = line
        if voffset_pix is not None and voffset_pix != 0:
            hdu.header['HISTORY'] = '--- Combination of AB spectra ---'
            hdu.header['HISTORY'] = "voffset_pix between A and B {}".format(
                voffset_pix)
        hdu.header.add_history('--- rinput.stored() (BEGIN) ---')
        for item in rinput.stored():
            value = getattr(rinput, item)
            cline = '{}: {}'.format(item, value)
            hdu.header.add_history(cline)
        hdu.header.add_history('--- rinput.stored() (END) ---')
        return result_img

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
