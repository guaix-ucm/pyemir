
import emirdrp.products as prod


def test_tags():
    m = prod.MasterRectWave().query_expr
    assert m.tags() == {'grism', 'filter'}
    m = prod.MasterIntensityFlat().query_expr
    assert m.tags() == {'filter'}