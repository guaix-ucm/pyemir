
import numpy

from ..utils import create_rot2d


def test_create_rot2d():
    import math
    angle =  math.pi / 2.0
    res = create_rot2d(angle)
    assert numpy.allclose(res, [[0, -1], [1, 0]])