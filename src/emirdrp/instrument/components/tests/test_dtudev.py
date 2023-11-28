import pytest

from ..dtu import DetectorTranslationUnit
from ..dtuaxis import DTUAxis
from emirdrp.instrument.tests.dtuheader import HEADERS


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuc1(hdr):
    dtu = DetectorTranslationUnit('DTU')
    assert len(dtu.children) == 3


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuc3(hdr):
    dtu = DetectorTranslationUnit('DTU')
    axis = dtu.get_device('DTU.X')
    assert isinstance(axis, DTUAxis)


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuc4(hdr):
    dtu = DetectorTranslationUnit('DTU')
    dtu.configure_with_header(hdr)
    val = dtu.get_property('DTU.X.coor')
    assert val == hdr['XDTU']


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuc2(hdr):
    dtu = DetectorTranslationUnit('DTU')
    dtu.configure_with_header(hdr)
    val = dtu.get_properties()
    assert sorted(val.keys()) == ['coor', 'coor_r', 'name', 'xdtu_r', 'ydtu_r', 'zdtu_r']


def test_dtuc5():
    dtu = DetectorTranslationUnit('DTU')
    info = {}
    key0 = 'DTU.X'
    val0 = -345.4
    info[key0] = {'coor': val0}
    dtu.configure(info)
    for node, conf in info.items():
        for key, val in conf.items():
            assert dtu.get_property(f'{node}.{key}') == val
