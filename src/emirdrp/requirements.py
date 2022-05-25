#
# Copyright 2008-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Typical requirements of recipes"""

from numina.core import Parameter, Requirement
from numina.core.requirements import ObservationResultRequirement
from numina.types.datatype import ListOfType

import emirdrp.products as prods


# Expose this from numina for convenience
EMIRObservationResultRequirement = ObservationResultRequirement


class EMIRConfigurationRequirement(Requirement):
    """The Recipe requires the configuration of EMIR."""
    def __init__(self):

        super(EMIRConfigurationRequirement,
              self).__init__(prods.EMIRConfigurationType,
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
              self).__init__(prods.MasterBadPixelMask,
                             'Master bad pixel mask',
                             optional=optional
                             )


class MasterBiasRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterBiasRequirement,
              self).__init__(prods.MasterBias,
                             'Master BIAS image',
                             optional=optional
                             )


class MasterDarkRequirement(Requirement):
    def __init__(self, optional=True):
        super(MasterDarkRequirement,
              self).__init__(prods.MasterDark,
                             'Master DARK image',
                             optional=optional
                             )


class MasterIntensityFlatFieldRequirement(Requirement):
    def __init__(self):
        super(MasterIntensityFlatFieldRequirement,
              self).__init__(prods.MasterIntensityFlat, 'Master intensity flatfield')


class MasterSpectralFlatFieldRequirement(Requirement):
    def __init__(self):
        super(MasterSpectralFlatFieldRequirement,
              self).__init__(prods.MasterSpectralFlat, 'Master spectral flatfield')


class RefinedBoundaryModelParamRequirement(Requirement):
    def __init__(self):
        super(RefinedBoundaryModelParamRequirement,
              self).__init__(prods.RefinedBoundaryModelParam, 'Refined boundary model')


class MasterRectWaveRequirement(Requirement):
    def __init__(self, optional=False):
        super(MasterRectWaveRequirement,
              self).__init__(prods.MasterRectWave, 'Master rectification+wavelength library',
                             optional=optional)


class RectWaveCoeffRequirement(Requirement):
    def __init__(self, optional=False):
        super(RectWaveCoeffRequirement,
              self).__init__(prods.RectWaveCoeff, 'Rectification+wavelength calibration coefficients',
                             optional=optional)


class ListOfRectWaveCoeffRequirement(Requirement):
    def __init__(self, optional=False):
        super(ListOfRectWaveCoeffRequirement,
              self).__init__(ListOfType(prods.RectWaveCoeff), 'List of Rectification+wavelength calibration coefficients',
                             optional=optional)


class MasterSkyRequirement(Requirement):
    def __init__(self, optional=False):
        super(MasterSkyRequirement,
              self).__init__(prods.MasterSky, 'Sky image for subtraction', optional=optional)


class SpectralSkyRequirement(Requirement):
    def __init__(self, optional=False):
        super(SpectralSkyRequirement,
              self).__init__(prods.SkySpectrum, 'Sky spectrum for subtraction', optional=optional)


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
