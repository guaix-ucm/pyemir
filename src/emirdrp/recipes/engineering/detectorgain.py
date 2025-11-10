#
# Copyright 2010-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Recipe for the reduction of gain calibration frames."""

import logging
import math

import numpy
import scipy.stats
from astropy.io import fits

from numina.core.requirements import ObservationResultRequirement
from numina.core import Parameter, DataFrame
from numina.exceptions import RecipeError
from numina.core import Result

from emirdrp.core.recipe import EmirRecipe
from emirdrp.instrument.channels import QUADRANTS
from emirdrp.instrument.channels import CHANNELS
from emirdrp.products import MasterGainMap, MasterRONMap

_logger = logging.getLogger("numina.recipes.emir")


class GainRecipe1(EmirRecipe):
    """Detector Gain Recipe.

    Recipe to calibrate the detector gain.
    """

    obresult = ObservationResultRequirement()
    region = Parameter(
        "channel",
        "Region used to compute: " "(full|quadrant|channel)",
        choices=["full", "quadrant", "channel"],
    )

    gain = Result(MasterGainMap(None, None, None))
    ron = Result(MasterRONMap(None, None))

    def region_method(self, reqs):
        mm = reqs["region"].tolower()
        if mm == "full":
            return (slice(0, 2048), slice(0, 2048))
        elif mm == "quadrant":
            return QUADRANTS
        elif mm == "channel":
            return CHANNELS
        else:
            raise ValueError

    def run(self, rinput):

        resets = []
        ramps = []

        for frame in rinput.obresult.frames:
            if frame.itype == "RESET":
                resets.append(frame.label)
                _logger.debug("%s is RESET", frame.label)
            elif frame.itype == "RAMP":
                ramps.append(frame.label)
                _logger.debug("%s is RAMP", frame.label)
            else:
                raise RecipeError("frame is neither a RAMP nor a RESET")

        channels = self.region_method(rinput)
        result_gain = numpy.zeros((len(channels),))
        result_ron = numpy.zeros_like(result_gain)

        counts = numpy.zeros((len(ramps), len(channels)))
        variance = numpy.zeros_like(counts)

        last_reset = resets[-1]
        _logger.debug("opening last reset image %s", last_reset)
        last_reset_data = fits.getdata(last_reset)

        for i, di in enumerate(ramps):
            with fits.open(di, mode="readonly") as fd:
                restdata = fd[0].data - last_reset_data
                for j, channel in enumerate(channels):
                    c = restdata[channel].mean()
                    _logger.debug("%f counts in channel", c)
                    counts[i, j] = c
                    v = restdata[channel].var(ddof=1)
                    _logger.debug("%f variance in channel", v)
                    variance[i, j] = v

        for j, _ in enumerate(channels):
            res = scipy.stats.linregress(counts[:, j], variance[:, j])
            slope, intercept, _r_value, _p_value, _std_err = res

            result_gain[j] = 1.0 / slope
            result_ron[j] = math.sqrt(intercept)
        cube = numpy.zeros((2, 2048, 2048))

        for gain, var, channel in zip(result_gain, result_ron, channels):
            cube[0][channel] = gain
            cube[1][channel] = var

        hdu = fits.PrimaryHDU(cube[0])
        hduvar = fits.ImageHDU(cube[1])
        hdulist = fits.HDUList([hdu, hduvar])

        gain = MasterGainMap(
            mean=result_gain, var=numpy.array([]), frame=DataFrame(hdulist)
        )
        ron = MasterRONMap(mean=result_ron, var=numpy.array([]))
        return self.create_result(gain=gain, ron=ron)
