
import numpy as np
from ..starcounts import SpagnaModel, RBModel, BSModel


def test_bsmodel():
    bsmodel = BSModel()
    assert np.allclose(bsmodel.integral_counts(18), 0.16328254083029348)
    assert np.allclose(bsmodel.differential_counts(18), 0.09013877275265489)


def test_rbmodel():
    rbmodel = RBModel()
    assert np.allclose(rbmodel.integral_counts(18), 0.7976500639286859)
    assert np.allclose(rbmodel.differential_counts(18), 0.11800000000000002)


def test_spagnamodel():
    sgmodel = SpagnaModel()
    assert np.allclose(sgmodel.integral_counts(18), 0.51772222)
    assert np.allclose(sgmodel.differential_counts(18), 5.4529320987654316e-05)
