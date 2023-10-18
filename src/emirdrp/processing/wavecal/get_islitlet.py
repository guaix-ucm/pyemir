#
# Copyright 2019-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
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
