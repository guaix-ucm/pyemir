#
# Copyright 2019 Universidad Complutense de Madrid
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

"""Return islitlet from Y-pixel coordinate in the rectified image"""

from emirdrp.core import EMIR_NBARS
from emirdrp.core import EMIR_NPIXPERSLIT_RECTIFIED


def get_islitlet(ipixel):
    """Return islilet (from 1 to EMIR_NBARS) from Y-pixel coordinate

    Parameters
    ----------
    ipixel : int
        Y-pixel coordinate in the rectified image.

    Returns
    -------
    islitlet : int
        Slitlet number (from 1 to EMIR_NBARS).

    """

    if ipixel < 1:
        raise ValueError('ipixel={} cannot be < 1'.format(ipixel))

    if ipixel > EMIR_NPIXPERSLIT_RECTIFIED * EMIR_NBARS:
        raise ValueError('ipixel={} cannot be > {}'.format(
            ipixel, EMIR_NPIXPERSLIT_RECTIFIED * EMIR_NBARS
        ))

    islitlet = int(ipixel / EMIR_NPIXPERSLIT_RECTIFIED) + 1

    if ipixel % EMIR_NPIXPERSLIT_RECTIFIED == 0:
        islitlet -= 1

    return islitlet
