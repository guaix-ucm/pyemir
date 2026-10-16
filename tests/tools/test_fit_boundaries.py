"""Checks of the boundaries of the slitlets of pyemir-fit_boundaries"""

import numpy

from emirdrp.tools.fit_boundaries import integrity_check


def boundaries(xdtu):
    """Boundaries of a slitlet in one exposure, with the DTU at xdtu"""
    return {
        "boundary_coef_lower": [1.0, 0.0],
        "boundary_coef_upper": [2.0, 0.0],
        "boundary_xmax_lower": 2000.0,
        "boundary_xmax_upper": 2000.0,
        "boundary_xmin_lower": 10.0,
        "boundary_xmin_upper": 10.0,
        "csu_bar_left": 100.0,
        "csu_bar_right": 110.0,
        "csu_bar_slit_center": 105.0,
        "csu_bar_slit_width": 10.0,
        "rotang": 0.0,
        "xdtu": xdtu,
        "xdtu_0": 0.0,
        "ydtu": 0.0,
        "ydtu_0": 0.0,
        "zdtu": 0.0,
        "zdtu_0": 0.0,
        "zzz_info1": "",
        "zzz_info2": "",
    }


def test_integrity_check():
    """The average and maximum differences of the DTU configurations of the exposures"""
    bounddict = {
        "meta_info": {"description": "slitlet boundaries from fits to continuum-lamp exposures"},
        "tags": {"grism": "J", "filter": "J"},
        "contents": {
            "slitlet01": {
                "2017-01-01T00:00:00.000": boundaries(1.0),
                "2017-01-01T00:10:00.000": boundaries(1.2),
            }
        },
    }
    average, maxdiff = integrity_check(bounddict, 1.0)
    assert numpy.isclose(average.xaxis.coor, 1.1)
    assert numpy.isclose(maxdiff.xaxis.coor, 0.2)
