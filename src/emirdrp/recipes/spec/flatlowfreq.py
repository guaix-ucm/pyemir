#
# Copyright 2019-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""
Spectroscopy mode, compute low-frequency flatfield
"""

import contextlib
import logging
import numpy as np
from scipy import ndimage
import uuid

import numina.array.combine as combine
from numina.array.display.ximplotxy import ximplotxy
from numina.array.display.pause_debugplot import pause_debugplot
from numina.array.wavecalib.fix_pix_borders import fix_pix_borders
from numina.array.wavecalib.apply_integer_offsets import apply_integer_offsets
from numina.core import Parameter
from numina.core import Result
from numina.frame.utils import copy_img
from numina.processing.combine import combine_imgs

from emirdrp.core.recipe import EmirRecipe
import emirdrp.datamodel
from emirdrp.instrument.csu_configuration import CsuConfiguration
from emirdrp.processing.wavecal.rectwv_coeff_from_mos_library \
    import rectwv_coeff_from_mos_library
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 import save_four_ds9
from emirdrp.processing.wavecal.rectwv_coeff_to_ds9 \
    import save_spectral_lines_ds9
from emirdrp.processing.wavecal.slitlet2d import Slitlet2D
import emirdrp.products as prods
import emirdrp.requirements as reqs
from emirdrp.tools.nscan_minmax_frontiers import nscan_minmax_frontiers

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS2
from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_MINIMUM_SLITLET_WIDTH_MM
from emirdrp.core import EMIR_MAXIMUM_SLITLET_WIDTH_MM


class SpecFlatLowFreq(EmirRecipe):
    """Process continuum exposures of continuum lamp (lamp ON-OFF)

    """

    logger = logging.getLogger(__name__)

    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_flat = reqs.MasterSpectralFlatFieldRequirement()
    rectwv_coeff = reqs.RectWaveCoeffRequirement(optional=True)
    master_rectwv = reqs.MasterRectWaveRequirement(optional=True)

    # note that 'sum' is not allowed as combination method
    method = Parameter(
        'sigmaclip',
        description='Combination method',
        choices=['mean', 'median', 'sigmaclip']
    )
    method_kwargs = Parameter(
        dict(),
        description='Arguments for combination method',
        optional=True
    )
    minimum_slitlet_width_mm = Parameter(
        float(EMIR_MINIMUM_SLITLET_WIDTH_MM),
        description='Minimum width (mm) for a valid slitlet',
        optional=True
    )
    maximum_slitlet_width_mm = Parameter(
        float(EMIR_MAXIMUM_SLITLET_WIDTH_MM),
        description='Maximum width (mm) for a valid slitlet',
        optional=True
    )
    global_integer_offset_x_pix = Parameter(
        0,
        description='Global offset (pixels) in wavelength direction (integer)',
        optional=True
    )
    global_integer_offset_y_pix = Parameter(
        0,
        description='Global offset (pixels) in spatial direction (integer)',
        optional=True
    )
    nwindow_x_median = Parameter(
        301,
        description='Window size to smooth the flatfield in the spectral '
                    'direction',
        optional=True
    )
    nwindow_y_median = Parameter(
        11,
        description='Window size to smooth the flatfield in the spatial '
                    'direction',
        optional=True
    )
    minimum_fraction = Parameter(
        0.01,
        description='Minimum useful flatfielding value (fraction of maximum)',
        optional=True
    )
    minimum_value_in_output = Parameter(
        0.01,
        description='Minimum value allowed in output: pixels below this value'
                    ' are set to 1.0',
        optional=True
    )
    maximum_value_in_output = Parameter(
        10.0,
        description='Maximum value allowed in output: pixels above this value'
                    ' are set to 1.0',
        optional=True
    )
    debugplot = Parameter(
        0,
        description='Debugger parameter',
        optional=True
    )

    reduced_flatlowfreq = Result(prods.MasterSpectralFlat)

    def run(self, rinput):
        self.logger.info('starting generation of flatlowfreq')

        self.logger.info('rectwv_coeff..........................: {}'.format(
            rinput.rectwv_coeff))
        self.logger.info('master_rectwv.........................: {}'.format(
            rinput.master_rectwv))
        self.logger.info('Minimum slitlet width (mm)............: {}'.format(
            rinput.minimum_slitlet_width_mm))
        self.logger.info('Maximum slitlet width (mm)............: {}'.format(
            rinput.maximum_slitlet_width_mm))
        self.logger.info('Global offset X direction (pixels)....: {}'.format(
            rinput.global_integer_offset_x_pix))
        self.logger.info('Global offset Y direction (pixels)....: {}'.format(
            rinput.global_integer_offset_y_pix))
        self.logger.info('nwindow_x_median......................: {}'.format(
            rinput.nwindow_x_median
        ))
        self.logger.info('nwindow_y_median......................: {}'.format(
            rinput.nwindow_y_median
        ))
        self.logger.info('Minimum fraction......................: {}'.format(
            rinput.minimum_fraction))
        self.logger.info('Minimum value in output...............: {}'.format(
            rinput.minimum_value_in_output))
        self.logger.info('Maximum value in output...............: {}'.format(
            rinput.maximum_value_in_output))

        # check rectification and wavelength calibration information
        if rinput.master_rectwv is None and rinput.rectwv_coeff is None:
            raise ValueError('No master_rectwv nor rectwv_coeff data have '
                             'been provided')
        elif rinput.master_rectwv is not None and \
                rinput.rectwv_coeff is not None:
            self.logger.warning('rectwv_coeff will be used instead of '
                                'master_rectwv')
        if rinput.rectwv_coeff is not None and \
                (rinput.global_integer_offset_x_pix != 0 or
                 rinput.global_integer_offset_y_pix != 0):
            raise ValueError('global_integer_offsets cannot be used '
                             'simultaneously with rectwv_coeff')

        # check headers to detect lamp status (on/off)
        list_lampincd = []
        for fname in rinput.obresult.frames:
            with fname.open() as f:
                list_lampincd.append(f[0].header['lampincd'])

        # Note: after installing the new Hawaii 2 RG detector, the
        # keyword LAMPINCD values are no longer 0 / 1 but FALSE / TRUE.
        # To avoid changing the code below, we replace
        # FALSE by 0 and TRUE by 1 when necessary
        contains_only_0_and_1 = all(item in [0, 1] for item in list_lampincd)
        if not contains_only_0_and_1:
            contains_only_TRUE_and_FALSE = all(item in ['TRUE', 'FALSE'] for item in list_lampincd)
            if not contains_only_TRUE_and_FALSE:
                self.logger.info(f'List of LAMPINCD values: {list_lampincd}')
                raise ValueError('Unexpected values of the keyword LAMPINCD')
            # replace (in place) FALSE by 0 and TRUE by 1 in this case
            list_lampincd = [1 if item == 'TRUE' else 0 for item in list_lampincd]

        # check number of images
        nimages = len(rinput.obresult.frames)
        n_on = list_lampincd.count(1)
        n_off = list_lampincd.count(0)
        self.logger.info('Number of images with lamp ON.........: {}'.format(
            n_on))
        self.logger.info('Number of images with lamp OFF........: {}'.format(
            n_off))
        self.logger.info('Total number of images................: {}'.format(
            nimages))
        if n_on == 0:
            raise ValueError('Insufficient number of images with lamp ON')
        if n_on + n_off != nimages:
            raise ValueError('Number of images does not match!')

        # check combination method
        if rinput.method_kwargs == {}:
            method_kwargs = None
        else:
            if rinput.method == 'sigmaclip':
                method_kwargs = rinput.method_kwargs
            else:
                raise ValueError('Unexpected method_kwargs={}'.format(
                    rinput.method_kwargs))

        # build object to proceed with bpm, bias, and dark (not flat)
        flow = self.init_filters(rinput)

        # available combination methods
        method = getattr(combine, rinput.method)

        # basic reduction of images with lamp ON or OFF
        lampmode = {0: 'off', 1: 'on'}
        reduced_image_on = None
        reduced_image_off = None
        for imode in lampmode.keys():
            self.logger.info('starting basic reduction of images with'
                             ' lamp {}'.format(lampmode[imode]))
            tmplist = [rinput.obresult.frames[i] for i, lampincd in
                       enumerate(list_lampincd) if lampincd == imode]
            if len(tmplist) > 0:
                with contextlib.ExitStack() as stack:
                    hduls = [stack.enter_context(fname.open()) for fname
                             in tmplist]
                    reduced_image = combine_imgs(
                        hduls,
                        method=method,
                        method_kwargs=method_kwargs,
                        errors=False,
                        prolog=None
                    )
                if imode == 0:
                    reduced_image_off = flow(reduced_image)
                    hdr = reduced_image_off[0].header
                    self.set_base_headers(hdr)
                    self.save_intermediate_img(reduced_image_off,
                                               'reduced_image_off.fits')
                elif imode == 1:
                    reduced_image_on = flow(reduced_image)
                    hdr = reduced_image_on[0].header
                    self.set_base_headers(hdr)
                    self.save_intermediate_img(reduced_image_on,
                                               'reduced_image_on.fits')
                else:
                    raise ValueError('Unexpected imode={}'.format(imode))

        # computation of ON-OFF
        header_on = reduced_image_on[0].header
        data_on = reduced_image_on[0].data.astype('float32')
        if n_off > 0:
            header_off = reduced_image_off[0].header
            data_off = reduced_image_off[0].data.astype('float32')
        else:
            header_off = None
            data_off = np.zeros_like(data_on)
        reduced_data = data_on - data_off

        # update reduced image header
        reduced_image = self.create_reduced_image(
            rinput,
            reduced_data,
            header_on,
            header_off,
            list_lampincd
        )

        # save intermediate image in work directory
        self.save_intermediate_img(reduced_image, 'reduced_image.fits')

        # define rectification and wavelength calibration coefficients
        if rinput.rectwv_coeff is None:
            rectwv_coeff = rectwv_coeff_from_mos_library(
                reduced_image,
                rinput.master_rectwv
            )
            # set global offsets
            rectwv_coeff.global_integer_offset_x_pix = \
                rinput.global_integer_offset_x_pix
            rectwv_coeff.global_integer_offset_y_pix = \
                rinput.global_integer_offset_y_pix
        else:
            rectwv_coeff = rinput.rectwv_coeff
        # save as JSON in work directory
        self.save_structured_as_json(rectwv_coeff, 'rectwv_coeff.json')
        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(rectwv_coeff)
            save_spectral_lines_ds9(rectwv_coeff)

        # apply global offsets (to both, the original and the cleaned version)
        image2d = apply_integer_offsets(
            image2d=reduced_data,
            offx=rectwv_coeff.global_integer_offset_x_pix,
            offy=rectwv_coeff.global_integer_offset_y_pix
        )

        # load CSU configuration
        mecs_header = emirdrp.datamodel.get_mecs_header(reduced_image)
        csu_conf = CsuConfiguration.define_from_header(
            mecs_header
        )
        # determine (pseudo) longslits
        dict_longslits = csu_conf.pseudo_longslits()

        # valid slitlet numbers
        list_valid_islitlets = list(range(1, EMIR_NBARS + 1))
        for idel in rectwv_coeff.missing_slitlets:
            self.logger.info('-> Removing slitlet (not defined): ' + str(idel))
            list_valid_islitlets.remove(idel)
        # filter out slitlets with widths outside valid range
        list_outside_valid_width = []
        for islitlet in list_valid_islitlets:
            slitwidth = csu_conf.csu_bar_slit_width(islitlet)
            if (slitwidth < rinput.minimum_slitlet_width_mm) or \
                    (slitwidth > rinput.maximum_slitlet_width_mm):
                list_outside_valid_width.append(islitlet)
                self.logger.info('-> Removing slitlet (width out of range): ' +
                                 str(islitlet))
        if len(list_outside_valid_width) > 0:
            for idel in list_outside_valid_width:
                list_valid_islitlets.remove(idel)

        # initialize rectified image
        image2d_flatfielded = np.zeros((EMIR_NAXIS2, EMIR_NAXIS1))

        # main loop
        # grism_name = rectwv_coeff.tags['grism']
        # filter_name = rectwv_coeff.tags['filter']
        cout = '0'
        debugplot = rinput.debugplot
        for islitlet in list(range(1, EMIR_NBARS + 1)):
            if islitlet in list_valid_islitlets:
                # define Slitlet2D object
                slt = Slitlet2D(islitlet=islitlet,
                                rectwv_coeff=rectwv_coeff,
                                debugplot=debugplot)
                if abs(slt.debugplot) > 10:
                    print(slt)

                # extract (distorted) slitlet from the initial image
                slitlet2d = slt.extract_slitlet2d(
                    image_2k2k=image2d,
                    subtitle='original image'
                )

                # rectify slitlet
                slitlet2d_rect = slt.rectify(
                    slitlet2d=slitlet2d,
                    resampling=2,
                    subtitle='original (cleaned) rectified'
                )
                naxis2_slitlet2d, naxis1_slitlet2d = slitlet2d_rect.shape

                if naxis1_slitlet2d != EMIR_NAXIS1:
                    print('naxis1_slitlet2d: ', naxis1_slitlet2d)
                    print('EMIR_NAXIS1.....: ', EMIR_NAXIS1)
                    raise ValueError("Unexpected naxis1_slitlet2d")

                # get useful slitlet region (use boundaries)
                spectrail = slt.list_spectrails[0]
                yy0 = slt.corr_yrect_a + \
                    slt.corr_yrect_b * spectrail(slt.x0_reference)
                ii1 = int(yy0 + 0.5) - slt.bb_ns1_orig
                spectrail = slt.list_spectrails[2]
                yy0 = slt.corr_yrect_a + \
                    slt.corr_yrect_b * spectrail(slt.x0_reference)
                ii2 = int(yy0 + 0.5) - slt.bb_ns1_orig

                # median spectrum
                sp_collapsed = np.median(slitlet2d_rect[ii1:(ii2 + 1), :],
                                         axis=0)

                # smooth median spectrum along the spectral direction
                sp_median = ndimage.median_filter(
                    sp_collapsed,
                    5,
                    mode='nearest'
                )

                ymax_spmedian = sp_median.max()
                y_threshold = ymax_spmedian * rinput.minimum_fraction
                lremove = np.where(sp_median < y_threshold)
                sp_median[lremove] = 0.0

                if abs(slt.debugplot) % 10 != 0:
                    xaxis1 = np.arange(1, naxis1_slitlet2d + 1)
                    title = 'Slitlet#' + str(islitlet) + ' (median spectrum)'
                    ax = ximplotxy(xaxis1, sp_collapsed,
                                   title=title,
                                   show=False,
                                   **{'label': 'collapsed spectrum'})
                    ax.plot(xaxis1, sp_median, label='fitted spectrum')
                    ax.plot([1, naxis1_slitlet2d], 2 * [y_threshold],
                            label='threshold')
                    # ax.plot(xknots, yknots, 'o', label='knots')
                    ax.legend()
                    ax.set_ylim(-0.05 * ymax_spmedian, 1.05 * ymax_spmedian)
                    pause_debugplot(slt.debugplot,
                                    pltshow=True, tight_layout=True)

                # generate rectified slitlet region filled with the
                # median spectrum
                slitlet2d_rect_spmedian = np.tile(sp_median,
                                                  (naxis2_slitlet2d, 1))
                if abs(slt.debugplot) % 10 != 0:
                    slt.ximshow_rectified(
                        slitlet2d_rect=slitlet2d_rect_spmedian,
                        subtitle='rectified, filled with median spectrum'
                    )

                # compute smooth surface
                # clipped region
                slitlet2d_rect_clipped = slitlet2d_rect_spmedian.copy()
                slitlet2d_rect_clipped[:(ii1-1), :] = 0.0
                slitlet2d_rect_clipped[(ii2+2):, :] = 0.0
                # unrectified clipped image
                slitlet2d_unrect_clipped = slt.rectify(
                    slitlet2d=slitlet2d_rect_clipped,
                    resampling=2,
                    inverse=True,
                    subtitle='unrectified, filled with median spectrum '
                             '(clipped)'
                )
                # normalize initial slitlet image (avoid division by zero)
                slitlet2d_norm_clipped = np.zeros_like(slitlet2d)
                for j in range(naxis1_slitlet2d):
                    for i in range(naxis2_slitlet2d):
                        den = slitlet2d_unrect_clipped[i, j]
                        if den == 0:
                            slitlet2d_norm_clipped[i, j] = 1.0
                        else:
                            slitlet2d_norm_clipped[i, j] = \
                                slitlet2d[i, j] / den
                # set to 1.0 one additional pixel at each side (since
                # 'den' above is small at the borders and generates wrong
                # bright pixels)
                slitlet2d_norm_clipped = fix_pix_borders(
                    image2d=slitlet2d_norm_clipped,
                    nreplace=1,
                    sought_value=1.0,
                    replacement_value=1.0
                )
                slitlet2d_norm_clipped = slitlet2d_norm_clipped.transpose()
                slitlet2d_norm_clipped = fix_pix_borders(
                    image2d=slitlet2d_norm_clipped,
                    nreplace=1,
                    sought_value=1.0,
                    replacement_value=1.0
                )
                slitlet2d_norm_clipped = slitlet2d_norm_clipped.transpose()
                slitlet2d_norm_smooth = ndimage.median_filter(
                    slitlet2d_norm_clipped,
                    size=(rinput.nwindow_y_median, rinput.nwindow_x_median),
                    mode='nearest'
                )

                if abs(slt.debugplot) % 10 != 0:
                    slt.ximshow_unrectified(
                        slitlet2d=slitlet2d_norm_clipped,
                        subtitle='unrectified, pixel-to-pixel (clipped)'
                    )
                    slt.ximshow_unrectified(
                        slitlet2d=slitlet2d_norm_smooth,
                        subtitle='unrectified, pixel-to-pixel (smoothed)'
                    )

                # ---

                # check for (pseudo) longslit with previous and next slitlet
                imin = dict_longslits[islitlet].imin()
                imax = dict_longslits[islitlet].imax()
                if islitlet > 1:
                    same_slitlet_below = (islitlet - 1) >= imin
                else:
                    same_slitlet_below = False
                if islitlet < EMIR_NBARS:
                    same_slitlet_above = (islitlet + 1) <= imax
                else:
                    same_slitlet_above = False

                for j in range(EMIR_NAXIS1):
                    xchannel = j + 1
                    y0_lower = slt.list_frontiers[0](xchannel)
                    y0_upper = slt.list_frontiers[1](xchannel)
                    n1, n2 = nscan_minmax_frontiers(y0_frontier_lower=y0_lower,
                                                    y0_frontier_upper=y0_upper,
                                                    resize=True)
                    # note that n1 and n2 are scans (ranging from 1 to NAXIS2)
                    nn1 = n1 - slt.bb_ns1_orig + 1
                    nn2 = n2 - slt.bb_ns1_orig + 1
                    image2d_flatfielded[(n1 - 1):n2, j] = \
                        slitlet2d_norm_smooth[(nn1 - 1):nn2, j]

                    # force to 1.0 region around frontiers
                    if not same_slitlet_below:
                        image2d_flatfielded[(n1 - 1):(n1 + 2), j] = 1
                    if not same_slitlet_above:
                        image2d_flatfielded[(n2 - 5):n2, j] = 1
                cout += '.'
            else:
                cout += 'i'

            if islitlet % 10 == 0:
                if cout != 'i':
                    cout = str(islitlet // 10)

            self.logger.info(cout)

        # restore global offsets
        image2d_flatfielded = apply_integer_offsets(
            image2d=image2d_flatfielded,
            offx=-rectwv_coeff.global_integer_offset_x_pix,
            offy=-rectwv_coeff.global_integer_offset_y_pix
        )

        # set pixels below minimum value to 1.0
        filtered = np.where(image2d_flatfielded <
                            rinput.minimum_value_in_output)
        image2d_flatfielded[filtered] = 1.0

        # set pixels above maximum value to 1.0
        filtered = np.where(image2d_flatfielded >
                            rinput.maximum_value_in_output)
        image2d_flatfielded[filtered] = 1.0

        # update image header
        reduced_flatlowfreq = self.create_reduced_image(
            rinput,
            image2d_flatfielded,
            header_on, header_off,
            list_lampincd
        )

        # ds9 region files (to be saved in the work directory)
        if self.intermediate_results:
            save_four_ds9(rectwv_coeff)
            save_spectral_lines_ds9(rectwv_coeff)

        # save results in results directory
        self.logger.info('end of flatlowfreq generation')
        result = self.create_result(
            reduced_flatlowfreq=reduced_flatlowfreq
        )
        return result

    def create_reduced_image(self, rinput, reduced_data,
                             header_on, header_off,
                             list_lampincd):
        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(fname.open()) for fname in
                     rinput.obresult.frames]
            # Copy header of first image
            # base_header = hduls[0][0].header.copy()
            result = copy_img(hduls[0])
            hdu = result[0]
            hdu.data = reduced_data
            self.set_base_headers(hdu.header)

            self.logger.debug('update result header')
            # update additional keywords
            hdu.header['UUID'] = str(uuid.uuid1())
            hdu.header['OBSMODE'] = 'flatlowfreq'
            hdu.header['TSUTC2'] = hduls[-1][0].header['TSUTC2']
            hdu.header['history'] = "Processed flatlowfreq"
            hdu.header['NUM-NCOM'] = (len(hduls), 'Number of combined frames')

            # update history
            dm = emirdrp.datamodel.EmirDataModel()
            for img, lampincd in zip(hduls, list_lampincd):
                imgid = dm.get_imgid(img)
                hdu.header['HISTORY'] = "Image '{}' has lampincd='{}'".format(
                    imgid, lampincd)
        hdu.header['HISTORY'] = "Processed flatlowfreq"
        hdu.header['HISTORY'] = '--- Reduction of images with lamp ON ---'
        for line in header_on['HISTORY']:
            hdu.header['HISTORY'] = line
        if header_off is not None:
            hdu.header['HISTORY'] = '--- Reduction of images with lamp OFF ---'
            for line in header_off['HISTORY']:
                hdu.header['HISTORY'] = line
        hdu.header.add_history('--- rinput.stored() (BEGIN) ---')
        for item in rinput.stored():
            value = getattr(rinput, item)
            cline = '{}: {}'.format(item, value)
            hdu.header.add_history(cline)
        hdu.header.add_history('--- rinput.stored() (END) ---')
        return result

    def set_base_headers(self, hdr):
        newhdr = super(SpecFlatLowFreq, self).set_base_headers(hdr)
        # Update EXP to 0
        newhdr['EXP'] = 0
        return newhdr
