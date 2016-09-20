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

"""
Spectroscopy mode, ABBA
"""


from numina.core import Product
from numina.core.requirements import ObservationResultRequirement
import numina.exceptions
import numina.core

from emirdrp.core import EmirRecipe
import emirdrp.products as prods
import emirdrp.processing.datamodel
from emirdrp.processing.combine import basic_processing


class BaseABBARecipe(EmirRecipe):
    """Process images in ABBA mode"""

    obresult = ObservationResultRequirement()
    spec_abba = Product(prods.DataFrameType)


    @classmethod
    def build_recipe_input(cls, obsres, dal, pipeline='default'):
        return cls.build_recipe_input_gtc(obsres, dal, pipeline=pipeline)

    @classmethod
    def build_recipe_input_gtc(cls, obsres, dal, pipeline='default'):
        cls.logger.debug('start recipe input builder')
        print(dir(obsres))
        stareImagesIds = obsres.stareSpectraIds
        cls.logger.debug('Stare Spectra images IDS: ', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['frame'])

        newOR = numina.core.ObservationResult()
        newOR.frames = stareImages
        newRI = cls.create_input(obresult=newOR)
        cls.logger.debug('end recipe input builder')
        return newRI

    def run(self, rinput):
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
        dm = emirdrp.processing.datamodel.EmirDataModel()
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
        hdu.header['EMIRUUID'] = uuid.uuid1().hex
        # Headers of last image
        hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
        result = fits.HDUList([hdu])
        return result
