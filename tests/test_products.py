import pytest

import emirdrp.products as prod


def test_tags():
    m = prod.MasterRectWave().query_expr
    assert m.tags() == {"grism", "filter"}
    m = prod.MasterIntensityFlat().query_expr
    assert m.tags() == {"filter"}


@pytest.mark.parametrize("cls", [prod.RectWaveCoeff, prod.MasterRectWave, prod.RefinedBoundaryModelParam])
def test_query_expr_grism_filter(cls):
    """The products are selected from the registry by grism and filter"""
    product = cls()
    assert product.query_expr.tags() == {"grism", "filter"}
    pred = product.query_expr.fill_placeholders(grism="J", filter="J")
    assert pred.eval(grism="J", filter="J") is True
    assert pred.eval(grism="H", filter="J") is False
