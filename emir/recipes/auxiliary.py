#
# Copyright 2011 Sergio Pascual
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

'''Auxiliary Recipes for EMIR'''

import logging
import time

import numpy
import pyfits

from numina import RecipeBase, Image, __version__
from numina import FITSHistoryHandler
from numina.recipes import Image, Keyword

__all__ = ['Recipe', 'BiasRecipe', 'DarkRecipe', 'FlatRecipe']

_logger = logging.getLogger('emir.recipes')

class Recipe(RecipeBase):
    '''Null recipe'''

    __requires__ = []
    __provides__ = []

    def __init__(self):
        super(Recipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    def run(self, block):

            return {'products': {} }



class BiasRecipe(RecipeBase):
    '''Process BIAS images and create MASTER_BIAS.'''

    __requires__ = []
    __provides__ = [Image('master_bias', comment='Master bias image')]

    def __init__(self):
        super(BiasRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    def run(self, rb):

        history_header = pyfits.Header()

        fh =  FITSHistoryHandler(history_header)
        fh.setLevel(logging.INFO)
        _logger.addHandler(fh)

    	_logger.info('starting bias reduction')

        images = rb.images

        cdata = []

        try:
            for image in images:
                hdulist = pyfits.open(image, memmap=True, mode='readonly')
                cdata.append(hdulist)

            _logger.info('stacking images')
            data = numpy.zeros(cdata[0]['PRIMARY'].data.shape, dtype='float32')
            for hdulist in cdata:
                data += hdulist['PRIMARY'].data

            data /= len(cdata)
            data += 2.0

            hdu = pyfits.PrimaryHDU(data, header=cdata[0]['PRIMARY'].header)
    
            # update hdu header with
            # reduction keywords
            hdr = hdu.header
            hdr.update('FILENAME', 'master_bias-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'BIAS', 'Image type')
            hdr.update('NUMTYP', 'MASTER_BIAS', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'BiasRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([hdu])

            _logger.info('bias reduction ended')

            # merge header with HISTORY log
            hdr.ascardlist().extend(history_header.ascardlist())    

            return {'products': {'master_bias': hdulist}}
        except OSError as error:
            return {'error' : {'exception': str(error)}}
        finally:
            _logger.removeHandler(fh)

            for hdulist in cdata:
                hdulist.close()

class DarkRecipe(RecipeBase):
    '''Process DARK images and provide MASTER_DARK. '''

    __requires__ = [Image('master_bias', comment='Master bias image')]
    __provides__ = [Image('master_dark', comment='Master dark image')]

    def __init__(self):
        super(DarkRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    def run(self, block):

        # HISTORY logger
        history_header = pyfits.Header()

        fh =  FITSHistoryHandler(history_header)
        fh.setLevel(logging.INFO)
        _logger.addHandler(fh)

    	_logger.info('starting dark reduction')

        try:
            _logger.info('subtracting bias %s', str(self.parameters['master_bias']))
            with pyfits.open(self.parameters['master_bias'], mode='readonly') as master_bias:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd['primary'].data
                        data -= master_bias['primary'].data
                

            _logger.info('stacking images from block %d', block.id)

            base = block.images[0]
           
            with pyfits.open(base, memmap=True) as fd:
                data = fd['PRIMARY'].data.copy()
                hdr = fd['PRIMARY'].header
           
            for image in block.images[1:]:
                with pyfits.open(image, memmap=True) as fd:
                    add_data = fd['primary'].data
                    data += add_data

            hdu = pyfits.PrimaryHDU(data, header=hdr)
    
            # update hdu header with
            # reduction keywords
            hdr = hdu.header
            hdr.update('FILENAME', 'master_dark-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'DARK', 'Image type')
            hdr.update('NUMTYP', 'MASTER_DARK', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'DarkRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([hdu])

            _logger.info('dark reduction ended')

            # merge final header with HISTORY log
            hdr.ascardlist().extend(history_header.ascardlist())    

            return {'products': {'master_dark': hdulist}}
        except OSError as error:
            return {'error' : {'exception': str(error)}}
        finally:
            _logger.removeHandler(fh)

class FlatRecipe(RecipeBase):
    '''Process FLAT images and provide MASTER_FLAT. '''

    __requires__ = [Image('master_bias'),
                    Image('master_dark')]
    __provides__ = [Image('master_flat')]

    def __init__(self):
        super(FlatRecipe, self).__init__(
                        author="Sergio Pascual <sergiopr@fis.ucm.es>",
                        version="0.1.0"
                )

    def run(self, block):

        # HISTORY logger
        history_header = pyfits.Header()

        fh =  FITSHistoryHandler(history_header)
        fh.setLevel(logging.INFO)
        _logger.addHandler(fh)

    	_logger.info('starting flat reduction')

        try:
            _logger.info('subtracting bias %s', str(self.parameters['master_bias']))
            with pyfits.open(self.parameters['master_bias'], mode='readonly') as master_bias:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd['primary'].data
                        data -= master_bias['primary'].data
                

            _logger.info('subtracting dark %s', str(self.parameters['master_dark']))
            with pyfits.open(self.parameters['master_dark'], mode='readonly') as master_dark:
                for image in block.images:
                    with pyfits.open(image, memmap=True) as fd:
                        data = fd['primary'].data
                        data -= master_dark['primary'].data


            _logger.info('stacking images from block %d', block.id)

            base = block.images[0]
           
            with pyfits.open(base, memmap=True) as fd:
                data = fd['PRIMARY'].data.copy()
                hdr = fd['PRIMARY'].header
           
            for image in block.images[1:]:
                with pyfits.open(image, memmap=True) as fd:
                    add_data = fd['primary'].data
                    data += add_data

            # Normalize flat to mean 1.0
            data[:] = 1.0

            hdu = pyfits.PrimaryHDU(data, header=hdr)

            # update hdu header with
            # reduction keywords
            hdr = hdu.header
            hdr.update('FILENAME', 'master_flat-%(block_id)d.fits' % self.environ)
            hdr.update('IMGTYP', 'FALT', 'Image type')
            hdr.update('NUMTYP', 'MASTER_FLAT', 'Data product type')
            hdr.update('NUMXVER', __version__, 'Numina package version')
            hdr.update('NUMRNAM', 'FlatRecipe', 'Numina recipe name')
            hdr.update('NUMRVER', self.__version__, 'Numina recipe version')

            hdulist = pyfits.HDUList([hdu])

            _logger.info('flat reduction ended')

            # merge final header with HISTORY log
            hdr.ascardlist().extend(history_header.ascardlist())    

            return {'products': {'master_flat': hdulist}}
        finally:
            _logger.removeHandler(fh)

