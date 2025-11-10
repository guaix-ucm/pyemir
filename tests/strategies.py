#
# Copyright 2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


import hypothesis.strategies as st

"""Hypothesis strategies for DTU header testing"""


@st.composite
def dtu_dict(draw):
    headerX = {
        "coor_f": draw(st.floats(min_value=0.8, max_value=1.2)),
        "coor": draw(st.floats(min_value=-210, max_value=-200)),
        "coor_0": draw(st.floats(min_value=-310, max_value=-290)),
    }
    headerY = {
        "coor_f": draw(st.floats(min_value=0.8, max_value=1.2)),
        "coor": draw(st.floats(min_value=-210, max_value=-200)),
        "coor_0": draw(st.floats(min_value=-30, max_value=18)),
    }
    headerZ = {
        "coor_f": draw(st.floats(min_value=0.8, max_value=1.2)),
        "coor": draw(st.floats(min_value=-480, max_value=-100)),
        "coor_0": draw(st.floats(min_value=-470, max_value=123)),
    }

    return {"X": headerX, "Y": headerY, "Z": headerZ}


@st.composite
def dtu_fits_dict(draw):
    mapping = draw(dtu_dict())
    header = {}
    for prop in ["coor", "coor_f", "coor_0"]:
        for axis in ["X", "Y", "Z"]:
            header[f"{axis}DTU"] = mapping[axis][prop]
    return header


@st.composite
def dtu_fits_header(draw):
    import astropy.io.fits as fits

    mapping = draw(dtu_fits_dict())
    header = fits.Header(mapping)
    return header
