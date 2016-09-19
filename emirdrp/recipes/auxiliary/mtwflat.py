#
# Copyright 2016 Universidad Complutense de Madrid
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

"""Twilight Flat Recipe for a list of frames in different filters"""


import numpy
import numpy.polynomial.polynomial as nppol
import astropy.io.fits as fits
from numina.array.combine import median
from numina.core import Product
from numina.core.types import ListOfType
from numina.core.requirements import ObservationResultRequirement

from emirdrp.processing.info import gather_info_frames
from emirdrp.core import EmirRecipe
from emirdrp.products import MasterIntensityFlat
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.processing.combine import basic_processing_with_combination_frames


class MultiTwilightFlatRecipe(EmirRecipe):
    """Create a list of twilight flats"""
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    twflatframes = Product(ListOfType(MasterIntensityFlat))

    def run(self, rinput):

        results = []
        self.logger.info('starting multiflat flat reduction')

        flow = self.init_filters(rinput)
        saturation = 45000.0

        iinfo = gather_info_frames(rinput.obresult.frames)
        image_groups = {}
        self.logger.info('group images by filter')
        for idx, info in enumerate(iinfo):
            filt = info['filter']
            if filt not in image_groups:
                self.logger.debug('new filter %s', filt)
                image_groups[filt] = []
            img = rinput.obresult.frames[idx]
            self.logger.debug('image %s in group %s', img, filt)
            image_groups[filt].append(img)

        for filt, frames in image_groups.items():
            self.logger.info('processing filter %s', filt)

            # Uncomment this line and comment the following
            # to revert to non-ramp
            # res = self.run_per_filter(frames, flow)
            res = self.run_per_filter_ramp(frames, saturation=saturation)

            results.append(res)

        self.logger.info('end multiflat flat reduction')
        result = self.create_result(twflatframes=results)

        return result

    def run_per_filter(self, frames, flow):

        errors = True
        self.logger.debug('using errors: %s', errors)
        hdulist = basic_processing_with_combination_frames(frames, flow,
                                                    method=median,
                                                    errors=errors)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        mm = hdulist[0].data.mean()
        self.logger.info('mean value of flat is: %f', mm)
        hdr['CCDMEAN'] = mm

        self.logger.debug('normalize image')
        hdulist[0].data /= mm
        if errors:
            self.logger.debug('normalize VAR extension')
            hdulist['variance'].data /= (mm * mm)

        return hdulist

    def run_per_filter_ramp(self, frames, saturation):
        errors = False
        nimages = len(frames)
        if nimages == 0:
            raise ValueError('len(images) == 0')

        median_frames = numpy.empty((nimages,))
        exptime_frames = []
        utc_frames = []
        # generate 3D cube
        bshape = (2048, 2048)
        flat_frames = numpy.empty((bshape[0], bshape[1], nimages))
        for idx, frame in enumerate(frames):
            image = frame.open()
            flat_frames[:, :, idx] = image['primary'].data
            exptime_frames.append(image[0].header['EXPTIME'])
            median_frames[idx] = numpy.median(image['primary'].data)
            utc_frames.append(image[0].header['UTC'])

            self.logger.debug(
                'image %d exptime %f median %f UTC %s',
                idx, exptime_frames[idx],
                median_frames[idx],
                utc_frames[-1]
            )

        # filter saturated images
        good_images = median_frames < saturation
        ngood_images = good_images.sum()
        if ngood_images < 2:
            self.logger.warning('We have only %d good images', ngood_images)
            # Reference image
            ref_image = frames[0].open()
            slope_scaled = numpy.ones(bshape) * exptime_frames[0]
        else:
            nsaturated = nimages - good_images.sum()
            if nsaturated > 0:
                self.logger.debug('we have %d images with median value over saturation (%f)',
                                  nsaturated , saturation)

            m = flat_frames[:,:, good_images]
            # Reshape array to obtain a 2D array
            m_r = m.reshape((bshape[0] * bshape[1], ngood_images))
            self.logger.debug('fitting slopes with mean-squares')
            ll = nppol.polyfit(median_frames[good_images], m_r.T, deg=1)
            slope = ll[1].reshape(bshape)
            base = ll[0].reshape(bshape)

            # First good frame
            index_of_first_good = numpy.nonzero(good_images)[0][0]
            ref_image = frames[index_of_first_good].open()
            slope_scaled = slope * exptime_frames[index_of_first_good]

        hdr = ref_image[0].header
        self.set_base_headers(hdr)

        hdu = fits.PrimaryHDU(data=slope_scaled, header=hdr)
        hdulist = fits.HDUList(hdu)

        mm = hdulist[0].data.mean()
        self.logger.info('mean value of flat is: %f', mm)
        hdr['CCDMEAN'] = mm

        self.logger.debug('normalize image')
        hdulist[0].data /= mm
        if errors:
            self.logger.debug('normalize VAR extension')
            hdulist['variance'].data /= (mm * mm)

        return hdulist
