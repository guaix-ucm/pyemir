#
# Copyright 2008-2026 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


import logging
import warnings

import numina.processing as proc
import numina.util.node
from numina.processing.getters import get_corrector_for, load_corrector

import emirdrp.core as c
from emirdrp.processing.flatfield import FlatFieldCorrector

_logger = logging.getLogger(__name__)

get_corrector_p = get_corrector_for("master_bpm", proc.BadPixelCorrector)
get_corrector_s = get_corrector_for("master_sky", proc.SkyCorrector)
get_corrector_d = get_corrector_for("master_dark", proc.DarkCorrector)
# The intensity and spectral flats, both in master_flat
get_corrector_f = get_corrector_for("master_flat", FlatFieldCorrector)
get_corrector_sf = get_corrector_for("master_flat", FlatFieldCorrector)


def get_corrector_b(rinput, meta, ins, datamodel):
    """Corrector of the bias, only with the read modes that require it"""
    iinfo = meta["obresult"]
    if iinfo:
        mode = iinfo[0]["readmode"]
        if mode.lower() in c.EMIR_BIAS_MODES:
            use_bias = True
            _logger.info("readmode is %s, bias required", mode)
        else:
            use_bias = False
            _logger.info("readmode is %s, bias not required", mode)
    else:
        use_bias = False

    if use_bias:
        return load_corrector(rinput, "master_bias", proc.BiasCorrector, datamodel)
    else:
        _logger.info("ignoring bias")
        return numina.util.node.IdNode()


def get_corrector_gen(rinput, datamodel, CorrectorClass, key):
    """The corrector of the calibration rinput.<key>

    .. deprecated:: 0.25
        Use numina.processing.getters.load_corrector.
    """
    warnings.warn(
        "get_corrector_gen is deprecated, it will be removed in 0.26, use numina.processing.getters.load_corrector",
        DeprecationWarning,
        stacklevel=2,
    )
    return load_corrector(rinput, key, CorrectorClass, datamodel)


def get_checker(rinput, meta, ins, datamodel):
    from emirdrp.processing.checkers import Checker

    return Checker()
