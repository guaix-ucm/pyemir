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


import logging

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
        logger = logging.getLogger('numina.recipes.emir')

        results = []
        logger.info('starting multiflat flat reduction')

        flow = self.init_filters(rinput)

        iinfo = gather_info_frames(rinput.obresult.frames)
        image_groups = {}
        logger.info('group images by filter')
        for idx, info in enumerate(iinfo):
            filt = info['filter']
            if filt not in image_groups:
                logger.debug('new filter %s', filt)
                image_groups[filt] = []
            img = rinput.obresult.frames[idx]
            logger.debug('image %s in group %s', img, filt)
            image_groups[filt].append(img)

        for filt, frames in image_groups.items():
            logger.info('processing filter %s', filt)

            res = self.run_per_filter(frames, flow)

            results.append(res)

        logger.info('end multiflat flat reduction')
        result = self.create_result(twflatframes=results)

        return result

    def run_per_filter(self, frames, flow):

        logger = logging.getLogger('numina.recipes.emir')

        errors = True
        logger.debug('using errors: %s', errors)
        hdulist = basic_processing_with_combination_frames(frames, flow,
                                                    method=median,
                                                    errors=errors)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)
        mm = hdulist[0].data.mean()
        logger.info('mean value of flat is: %f', mm)
        hdr['CCDMEAN'] = mm

        logger.debug('normalize image')
        hdulist[0].data /= mm
        if errors:
            logger.debug('normalize VAR extension')
            hdulist['variance'].data /= (mm * mm)

        return hdulist
