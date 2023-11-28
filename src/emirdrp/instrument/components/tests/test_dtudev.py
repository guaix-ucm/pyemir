
from hypothesis import given
import hypothesis.strategies as st

from emirdrp.instrument.tests.dtuheader import dtu_dict, dtu_fits_header

from ..dtu import DetectorTranslationUnit
from ..dtuaxis import DTUAxis


def test_dtuc1():
    dtu = DetectorTranslationUnit('DTU')
    assert len(dtu.children) == 3


@given(st.sampled_from(['X', 'Y', 'Z']))
def test_dtuc2(name):
    dtu = DetectorTranslationUnit('DTU')
    axis = dtu.get_device(f'DTU.{name}')
    assert isinstance(axis, DTUAxis)


@given(dtu_fits_header(), st.sampled_from(['X', 'Y', 'Z']))
def test_dtuc3(hdr, name):
    dtu = DetectorTranslationUnit('DTU')
    dtu.configure_with_header(hdr)
    val = dtu.get_property(f'DTU.{name}.coor')
    assert val == hdr[f'{name}DTU']


@given(dtu_fits_header())
def test_dtuc4(hdr):
    dtu = DetectorTranslationUnit('DTU')
    dtu.configure_with_header(hdr)
    val = dtu.get_properties()
    assert sorted(val.keys()) == ['coor', 'coor_r', 'name', 'xdtu_r', 'ydtu_r', 'zdtu_r']


@given(dtu_dict())
def test_dtuc5(info):
    dtu = DetectorTranslationUnit('DTU')
    dtu.configure(info)
    for node, conf in info.items():
        for key, val in conf.items():
            assert dtu.get_property(f'{node}.{key}') == val
