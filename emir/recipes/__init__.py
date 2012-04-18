#
# Copyright 2008-2012 Universidad Complutense de Madrid
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

'''Recipes for EMIR Observing Modes.'''

from numina.pipeline import Pipeline

from .auxiliary import BiasRecipe, DarkRecipe, IntensityFlatRecipe
from .image import StareImageRecipe, DitheredImageRecipe, MicroditheredImageRecipe
from .image import NBImageRecipe, MosaicRecipe

_equiv_class = {
    'bias_image': BiasRecipe,
    'dark_image': DarkRecipe,
    'intensity_flatfield': IntensityFlatRecipe,
    'msm_spectral_flatfield': 'auxiliary:Recipe',
    'slit_transmission_calibration': 'auxiliary:Recipe',
    'wavelength_calibration': 'auxiliary:Recipe',
    'ts_rough_focus': 'auxiliary:Recipe',
    'ts_fine_focus': 'auxiliary:Recipe',
    'emir_focus_control': 'auxiliary:Recipe',
    'image_setup': 'auxiliary:Recipe',
    'mos_and_longslit_setup': 'auxiliary:Recipe',
    'target_acquisition': 'auxiliary:Recipe',
    'mask_imaging': 'auxiliary:Recipe',
    'msm_and_lsm_check': 'auxiliary:Recipe',
    'stare_image': StareImageRecipe,
    'nb_image': NBImageRecipe,
    'dithered_image': DitheredImageRecipe,
    'microdithered_image': MicroditheredImageRecipe,
    'mosaiced_image': MosaicRecipe,
    'stare_spectra': 'mos:Recipe',
    'dn_spectra': 'mos:Recipe',
    'offset_spectra': 'mos:Recipe',
    'raster_spectra': 'ls:Recipe',
}

class EmirPipeline(Pipeline):
    def __init__(self, version):
        super(EmirPipeline, self).__init__(name='emir', 
                version=version,
                recipes=_equiv_class)

