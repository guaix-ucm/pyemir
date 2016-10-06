#
# Copyright 2011-2015 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published byXS
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

"""MOS Recipes for EMIR"""

from numina.core import Parameter, Requirement
from numina.core import Product

from emirdrp.core import EmirRecipe
from numina.core.requirements import ObservationResultRequirement
from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
import emirdrp.products as prods


class DNSpectraRecipe(EmirRecipe):
    '''
    Observing mode:
        Dithered/Nodded spectra along the slit
    '''

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_spectral_ff = Requirement(prods.MasterSpectralFlat, 'Master spectral flatfield')
    st_calibration = Requirement(prods.SlitTransmissionCalibration, 'Slit tranmision calibration')
    w_calibration = Requirement(prods.WavelengthCalibration, 'Wavelength calibration')
    lines = Parameter('lines', None, 'List of x-lambda pairs of line coordinates')

    spectra = Product(prods.Spectra)
    catalog = Product(prods.LinesCatalog)

    def run(self, rinput):
        return self.create_result(
            spectra=prods.Spectra(),
            catalog=prods.LinesCatalog()
        )


class OffsetSpectraRecipe(EmirRecipe):
    """
    Observing mode:
        Offset spectra beyond the slit
    """

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_spectral_ff = Requirement(prods.MasterSpectralFlat, 'Master spectral flatfield')
    st_calibration = Requirement(prods.SlitTransmissionCalibration, 'Slit tranmision calibration')
    w_calibration = Requirement(prods.WavelengthCalibration, 'Wavelength calibration')
    lines = Parameter('lines', None, 'List of x-lambda pairs of line coordinates')

    spectra = Product(prods.Spectra)
    catalog = Product(prods.LinesCatalog)

    def run(self, rinput):
        return self.create_result(
            spectra=prods.Spectra(),
            catalog=prods.LinesCatalog()
        )


class RasterSpectraRecipe(EmirRecipe):
    """
    Observing mode:
        Raster spectra
    """

    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_spectral_ff = Requirement(prods.MasterSpectralFlat, 'Master spectral flatfield')
    st_calibration = Requirement(prods.SlitTransmissionCalibration, 'Slit tranmision calibration')
    w_calibration = Requirement(prods.WavelengthCalibration, 'Wavelength calibration')

    cube = Product(prods.DataCube)

    def run(self, rinput):
        return self.create_result(cube=prods.DataCube())
