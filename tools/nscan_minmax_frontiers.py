from __future__ import division
from __future__ import print_function

from emirdrp.core import EMIR_NAXIS2


def nscan_minmax_frontiers(y0_frontier_lower, y0_frontier_upper,
                           resize=False):
    """Compute valid scan range for provided y0_frontier values.

    Parameters
    ----------
    y0_frontier_lower : float
        Ordinate of the lower frontier.
    y0_frontier_upper : float
        Ordinate of the upper frontier.
    resize : bool
        If True, when the limits are beyond the expected values
        [1,EMIR_NAXIS2], the values are truncated.

    Returns
    -------
    nscan_min : int
        Minimum useful scan for the image.
    nscan_max : int
        Maximum useful scan for the image.

    """

    fraction_pixel = y0_frontier_lower - int(y0_frontier_lower)
    if fraction_pixel > 0.0:
        nscan_min = int(y0_frontier_lower) + 1
    else:
        nscan_min = int(y0_frontier_lower)
    if nscan_min < 1:
        if resize:
            nscan_min = 1
        else:
            raise ValueError("nscan_min=" + str(nscan_min) + " is < 1")

    fraction_pixel = y0_frontier_upper - int(y0_frontier_upper)
    if fraction_pixel > 0.0:
        nscan_max = int(y0_frontier_upper)
    else:
        nscan_max = int(y0_frontier_upper) - 1
    if nscan_max > EMIR_NAXIS2:
        if resize:
            nscan_max = EMIR_NAXIS2
        else:
            raise ValueError("nscan_max=" + str(nscan_max) +
                             " is > NAXIS2_EMIR=" + str(EMIR_NAXIS2))

    return nscan_min, nscan_max
