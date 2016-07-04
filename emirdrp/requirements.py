#
# Copyright 2008-2014 Universidad Complutense de Madrid
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

"""Typical requirements of recipes"""

from numina.core import Parameter, Requirement

from .products import MasterBias, MasterDark, MasterBadPixelMask
from .products import MasterIntensityFlat
from .products import MasterSpectralFlat
from .products import MasterSky
from .products import EMIRConfigurationType
# from .products import SourcesCatalog


class EMIRConfigurationRequirement(Requirement):
    """The Recipe requires the configuration of EMIR."""
    def __init__(self):

        super(EMIRConfigurationRequirement,
              self).__init__(EMIRConfigurationType,
                             "EMIR Configuration"
                             )


    def __repr__(self):
        sclass = type(self).__name__
        return "%s(dest=%r, description='%s')" % (sclass,
                                                  self.dest,
                                                  self.description
                                                  )


class MasterBadPixelMaskRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterBadPixelMaskRequirement,
              self).__init__(MasterBadPixelMask,
                             'Master bad pixel mask',
                             optional=optional
                             )


class MasterBiasRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterBiasRequirement,
              self).__init__(MasterBias,
                             'Master BIAS image',
                             optional=optional
                             )


class MasterDarkRequirement(Requirement):
    def __init__(self):
        super(MasterDarkRequirement,
              self).__init__(MasterDark, 'Master DARK image')


class MasterIntensityFlatFieldRequirement(Requirement):
    def __init__(self):
        super(MasterIntensityFlatFieldRequirement,
              self).__init__(MasterIntensityFlat, 'Master intensity flatfield')


class MasterSpectralFlatFieldRequirement(Requirement):
    def __init__(self):
        super(MasterSpectralFlatFieldRequirement,
              self).__init__(MasterSpectralFlat, 'Master spectral flatfield')


class MasterSkyRequirement(Requirement):
    def __init__(self, optional=False):
        super(MasterSkyRequirement,
              self).__init__(MasterSky, 'Sky image for subtraction', optional=optional)


class Extinction_Requirement(Parameter):
    def __init__(self):
        super(Extinction_Requirement, self).__init__(
            0.0, 'Mean atmospheric extinction'
            )


class SkyImageSepTime_Requirement(Parameter):
    def __init__(self):
        super(SkyImageSepTime_Requirement, self).__init__(
            10.0,
            'Maximum time interval between target and sky images [minutes]'
            )


class Catalog_Requirement(Parameter):
    def __init__(self, optional=True):
        super(Catalog_Requirement, self).__init__(
            [], 'List of x, y coordinates to measure FWHM',
            optional=optional
            )


class Offsets_Requirement(Parameter):
    def __init__(self, optional=True):
        super(Offsets_Requirement, self).__init__(
            [], 'List of pairs of offsets',
            optional=optional
            )
