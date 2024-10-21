#
# Copyright 2011-2024 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


"""
Stare Image mode of EMIR
"""

from copy import deepcopy
import datetime
import logging
import sys

from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
from astropy.wcs import WCS
import numpy as np
from reproject import reproject_interp, reproject_adaptive, reproject_exact
from scipy.ndimage import median_filter

from numina.array import combine
from numina.core import Result, Parameter
from numina.core.query import Ignore
from numina.core.recipes import timeit
from numina.processing.combine import basic_processing_with_combination
from numina.util.context import manage_fits

from emirdrp.core.recipe import EmirRecipe
import emirdrp.core.extra as extra
import emirdrp.requirements as reqs
import emirdrp.products as prods
import emirdrp.processing.combine as comb
import emirdrp.decorators


if sys.version_info[:2] <= (3, 10):
    datetime.UTC = datetime.timezone.utc


_logger = logging.getLogger(__name__)


class StareImageRecipe2(EmirRecipe):
    """Process images in Stare Image Mode"""

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    reprojection_method = Parameter(
        'interp',
        description='Astrometric reprojection method',
        choices=['interp', 'adaptive', 'exact', 'none']
    )

    reduced_image = Result(prods.ProcessedImage)

    def run(self, rinput):
        self.logger.info('starting stare image reduction (offline)')

        # determine which EMIR detector we are using
        insconf = rinput.obresult.configuration
        detector_channels = insconf.get_device("detector").get_property("channels")
        self.logger.info(f"Detector channels: {detector_channels}")

        frames = rinput.obresult.frames

        with manage_fits(frames) as list_of:
            c_img = comb.combine_images(
                list_of, method=combine.mean, errors=False)

        flow = self.init_filters(rinput)

        # Correct Bias if needed
        # Correct Dark if needed
        # Correct FF

        processed_img = flow(c_img)

        hdr = processed_img[0].header
        self.set_base_headers(hdr)

        if rinput.master_bpm is not None:
            self.logger.debug('using BPM from inputs')
            hdul_bpm = rinput.master_bpm.open()
            hdu_bpm = extra.generate_bpm_hdu(hdul_bpm[0])
        else:
            self.logger.debug('using empty BPM')
            hdu_bpm = extra.generate_empty_bpm_hdu(processed_img[0])

        # remove ERROR extension
        if 'ERROR' in processed_img:
            self.logger.debug('... removing ERROR extension')
            processed_img.pop('ERROR')

        # image reprojection: note that the reprojection functions preserve
        # the flux only when the image is in surface brightness units
        convert_to_surface_brightness = True
        reprojection_method = rinput.reprojection_method
        if reprojection_method != 'none':
            if reprojection_method == 'interp':
                reproject_function = reproject_interp
            elif reprojection_method == 'adaptive':
                reproject_function = reproject_adaptive
            elif reprojection_method == 'exact':
                reproject_function = reproject_exact
            else:
                raise ValueError(
                    f'Unexpected astrometric_reprojection value: {reprojection_method}')
            # reproject data
            self.logger.debug('starting image reprojection')
            hdr_original = deepcopy(hdr)
            wcs_original = WCS(hdr_original)
            # remove PV2_2, PV2_3,... PV2_5
            for item in wcs_original.wcs.get_pv():
                if item[1] == 1:
                    self.logger.debug(
                        f'... preserving keyword PV{item[0]}_{item[1]}={item[2]}')
                elif item[1] > 1:
                    self.logger.debug(
                        f'... removing   keyword PV{item[0]}_{item[1]}={item[2]}')
                    del hdr[f'PV{item[0]}_{item[1]}']
                else:
                    raise ValueError('Unexpected PVi_j value')
            wcs_final = WCS(hdr)
            naxis2, naxis1 = processed_img[0].data.shape
            if convert_to_surface_brightness:
                # solid angle subtended by every pixel
                self.logger.debug('... computing solid angle of every pixel in the original WCS')
                pixel_solid_angle_original = pixel_solid_angle_arcsec2(
                    wcs=wcs_original,
                    naxis1=naxis1,
                    naxis2=naxis2,
                    kernel_size=(11, 11)
                )
                surface_brightness_data = processed_img[0].data / \
                    pixel_solid_angle_original
            else:
                surface_brightness_data = processed_img[0].data
            # reprojection itself
            self.logger.debug(
                f'... reprojecting surface brightness using reproject_{reprojection_method}')
            data_final, footprint = reproject_function(
                input_data=(surface_brightness_data, wcs_original),
                output_projection=wcs_final,
                shape_out=processed_img[0].data.shape,
            )
            if convert_to_surface_brightness:
                self.logger.debug('... computing solid angle of every pixel in the final WCS')
                pixel_solid_angle_final = pixel_solid_angle_arcsec2(
                    wcs=wcs_final,
                    naxis1=naxis1,
                    naxis2=naxis2,
                    kernel_size=(11, 11)
                )
                data_final *= pixel_solid_angle_final
            # avoid undefined values: note that using footprint < 1.0 does not work
            # properly with reproject_exact(), which gives a footprint with values around 1.0
            # but with a non-negligible dispersion below 0.01 (for safety, here we use
            # a threshold of 0.98)
            minimum_footprint = 0.98
            # avoid undefined values
            data_final[footprint < minimum_footprint] = 0
            processed_img[0].data = data_final
            mask_footprint = (footprint < minimum_footprint).astype('uint8')
            # reproject mask
            self.logger.debug(
                f'... reprojecting mask using reproject_{reprojection_method}')
            mask_reprojected, footprint = reproject_function(
                input_data=(hdu_bpm.data, wcs_original),
                output_projection=wcs_final,
                shape_out=hdu_bpm.data.shape,
            )
            # avoid undefined values
            mask_reprojected[footprint < minimum_footprint] = 0
            self.logger.debug(
                '... merging existing mask with footprint from reproject')
            # there is no need to recompute mask_footprint (is the one computed for data)
            hdu_bpm.data = (mask_reprojected > 0).astype(
                'uint8') + mask_footprint
            # generate image with individual channel distortion
            if detector_channels == 'H2RG_FULL':
                image_channels = np.zeros((2048, 2048))
                for j_channel in range(32):
                    j1 = j_channel * 64
                    j2 = j1 + 64
                    image_channels[:, j1:j2] = j_channel + 1
                channels_reprojected, footprint = reproject_function(
                    input_data=(image_channels, wcs_original),
                    output_projection=wcs_final,
                    shape_out=image_channels.shape,
                )
                channels_reprojected[footprint < minimum_footprint] = 0.0
                channels_reprojected_int = np.round(channels_reprojected).astype('uint8')
                header = fits.Header()
                header['EXTNAME'] = 'ICHANNEL'
                hdu_channels = fits.ImageHDU(channels_reprojected_int, header=header)
            else:
                hdu_channels = None

            # update header
            hdr['history'] = f'Reprojection using reproject_{reprojection_method}()'
            hdr['history'] = f'Reprojection time {datetime.datetime.now(datetime.UTC).isoformat()}'
        else:
            if detector_channels == 'H2RG_FULL':
                image_channels = np.zeros((2048, 2048), dtype='uint8')
                for j_channel in range(32):
                    j1 = j_channel * 64
                    j2 = j1 + 64
                    image_channels[:, j1:j2] = j_channel + 1
                header = fits.Header()
                header['EXTNAME'] = 'ICHANNEL'
                hdu_channels = fits.ImageHDU(image_channels, header=header)
            else:
                hdu_channels = None

        # Append the BPM to the result
        self.logger.debug('append BPM')
        processed_img.append(hdu_bpm)

        # Append the reprojected layout of the detector channels (for the H2RG detector)
        if hdu_channels is not None:
            self.logger.debug('append channels reprojection')
            processed_img.append(hdu_channels)

        self.logger.info('end stare image (off) reduction')
        result = self.create_result(reduced_image=processed_img)

        return result

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(StareImageRecipe2, self).set_base_headers(hdr)
        # Set EXP to 0
        hdr['EXP'] = 0
        return hdr


class StareImageBaseRecipe(EmirRecipe):
    """Process images in Stare Image Mode"""

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement(optional=True)

    frame = Result(prods.ProcessedImage)

    def __init__(self, *args, **kwargs):
        super(StareImageBaseRecipe, self).__init__(*args, **kwargs)
        if False:
            self.query_options['master_sky'] = Ignore()

    @emirdrp.decorators.loginfo
    @timeit
    def run(self, rinput):
        self.logger.info('starting stare image reduction')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(
            rinput,
            flow,
            method=combine.median
        )
        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        if rinput.master_bpm:
            hdul_bpm = rinput.master_bpm.open()
            hdu_bpm = extra.generate_bpm_hdu(hdul_bpm[0])
        else:
            hdu_bpm = extra.generate_empty_bpm_hdu(hdulist[0])

        # Append the BPM to the result
        hdulist.append(hdu_bpm)
        self.logger.info('end stare image reduction')
        result = self.create_result(frame=hdulist)

        return result

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(StareImageBaseRecipe, self).set_base_headers(hdr)
        # Update EXP to 0
        hdr['EXP'] = 0
        return hdr


def pixel_solid_angle_arcsec2(wcs, naxis1, naxis2, kernel_size=None):
    """Compute the solid angle (arcsec**2) of every pixel."""

    # X, Y coordinates (2D image array, following the FITS criterium)
    # corresponding to the four corners of all the image pixels
    borders_x = np.arange(naxis1 + 1) + 0.5
    borders_y = np.arange(naxis2 + 1) + 0.5
    meshgrid = np.meshgrid(borders_x, borders_y)
    ix_array = meshgrid[0].flatten()
    iy_array = meshgrid[1].flatten()

    # spherical coordinates of the four corners of all the image pixels
    result_spherical = SkyCoord.from_pixel(
        xp=ix_array,
        yp=iy_array,
        wcs=wcs,
        origin=1,
        mode='all'
    )

    # cartesian coordinates of the four corners of all the image pixels
    x = result_spherical.cartesian.x.value.reshape(naxis2 + 1, naxis1 + 1)
    y = result_spherical.cartesian.y.value.reshape(naxis2 + 1, naxis1 + 1)
    z = result_spherical.cartesian.z.value.reshape(naxis2 + 1, naxis1 + 1)

    # dot product of consecutive points along NAXIS1
    dot_product_naxis1 = x[:, :-1] * x[:, 1:] + \
        y[:, :-1] * y[:, 1:] + z[:, :-1] * z[:, 1:]
    # distance (arcsec) between consecutive points along NAXIS1
    result_naxis1 = np.arccos(dot_product_naxis1) * 180 / np.pi * 3600
    # average distances corresponding to the upper and lower sides of each pixel
    pixel_size_naxis1 = (result_naxis1[:-1, :] + result_naxis1[1:, :]) / 2

    # dot product of consecutive points along NAXIS2
    dot_product_naxis2 = x[:-1, :] * x[1:, :] + \
        y[:-1, :] * y[1:, :] + z[:-1, :] * z[1:, :]
    # distance (arcsec) between consecutive points along NAXIS2
    result_naxis2 = np.arccos(dot_product_naxis2) * 180 / np.pi * 3600
    # averange distances corresponding to the left and right sides of each pixel
    pixel_size_naxis2 = (result_naxis2[:, :-1] + result_naxis2[:, 1:]) / 2

    # pixel size (arcsec**2)
    result = pixel_size_naxis1 * pixel_size_naxis2

    # smooth result
    if kernel_size is not None:
        result = median_filter(result, size=kernel_size, mode='nearest')

    return result
