#
# Copyright 2008-2011 Sergio Pascual
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

'''Recipes for Emir Observing Modes.'''

__all__ = ['find_recipe']

class EmirRecipeMixin:
    instrument = ['emir']

# equivalence
_equiv = {
    'bias_image': 'auxiliary:Recipe',
    'dark_current_image': 'auxiliary:Recipe',
    'intensity_flatfield': 'auxiliary:Recipe',
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
    'stare_image': 'image:Recipe',
    'nb_image': 'image:Recipe',
    'dithered_image':'image:Recipe',
    'microdithered_image':'image:Recipe',
    'mosaiced_image': 'image:Recipe',
    'stare_spectra': 'mos:Recipe',
    'dn_spectra': 'mos:Recipe',
    'offset_spectra': 'mos:Recipe',
    'raster_spectra': 'ls:Recipe',
    
}

def find_recipe(mode):
    return _equiv[mode]
