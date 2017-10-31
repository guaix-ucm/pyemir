#
# Copyright 2011-2016 Universidad Complutense de Madrid
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
Stare Image mode of EMIR
"""

import numpy
import astropy.io.fits as fits
from numina.array import combine
from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
from numina.core.query import Ignore

from emirdrp.core import EmirRecipe
from emirdrp.products import DataFrameType
import emirdrp.requirements as reqs
from emirdrp.processing.combine import basic_processing_with_combination
import emirdrp.decorators


class StareImageBaseRecipe(EmirRecipe):
    """Process images in Stare Image Mode"""

    obresult = ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_flat = reqs.MasterIntensityFlatFieldRequirement()
    master_sky = reqs.MasterSkyRequirement(optional=True)

    frame = Product(DataFrameType)

    def __init__(self, *args, **kwargs):
        super(StareImageBaseRecipe, self).__init__(*args, **kwargs)
        if False:
            self.query_options['master_sky'] = Ignore()

    @emirdrp.decorators.loginfo
    @emirdrp.decorators.timeit
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
            hdu_bpm = generate_bpm_hdu(hdul_bpm[0])
        else:
            hdu_bpm = generate_empty_bpm_hdu(hdulist[0])

        # Append the BPM to the result
        hdulist.append(hdu_bpm)
        self.logger.info('end stare image reduction')
        result = self.create_result(frame=hdulist)

        return result

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(StareImageBaseRecipe, self).set_base_headers(hdr)
        # Update SEC to 0
        hdr['SEC'] = 0
        return hdr


def generate_bpm_hdu(hdu):
    return generate_bpm_data(hdu.data, hdu.header)


def generate_empty_bpm_hdu(hdu):
    data = numpy.zeros_like(hdu.data, dtype='uint8')
    return generate_bpm_data(data)


def generate_bpm_data(data, header=None):
    hdu_bpm = fits.ImageHDU(
        data,
        header=header
    )

    hdu_bpm.header['EXTNAME'] = 'BPM'
    return hdu_bpm
