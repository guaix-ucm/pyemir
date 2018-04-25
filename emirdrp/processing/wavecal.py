#
# Copyright 2018 Universidad Complutense de Madrid
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


"""Wavelength calibration functionality"""

import logging


def rect_wavecal(reduced_image, master_rectwv):
    """Apply rectification and wavelength calibration.

    Parameters
    ----------
    reduced_image : HDUList object
        Image with preliminary basic reduction: bpm,

    Returns
    -------
    TBD

    """

    logger = logging.getLogger(__name__)

    print(type(reduced_image))
    print(type(master_rectwv))

    logger.info("something")
    logger.debug("something else")
    logger.warning("warning %s", "mensaje")
