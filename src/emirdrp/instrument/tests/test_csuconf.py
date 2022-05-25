
import pkgutil

import pytest

from ..csuconf import CSUConf
from ..csuconf import CSUBarModel, CSUBarModelL, CSUBarModelR
from ..csuconf import PhysicalBar, PhysicalBarL, PhysicalBarR
from ..csuconf import LogicalSlit, TargetType
from ..csuconf import create_bar_models, read_csu_from_header
from ..csuconf import merge_slits
from ..csuconf import EMIR_NBARS


import numpy
try:
    import StringIO as S
except ImportError:
    import io as S


def create_test_header0():
    hdr = {}
    for i in range(55):
        hdr["CSUP{}".format(i + 1)] = -100
    for i in range(55, 110):
        hdr["CSUP{}".format(i + 1)] = -100
    return hdr


def create_test_header1():
    hdr = create_test_header0()

    hdr["SLIFL12"] = 2
    hdr["SLIFL13"] = 2
    hdr["SLIFL33"] = 2
    hdr["SLIFL44"] = 2
    hdr["SLIFL45"] = 2

    hdr["XRSLI12"] = 200.1
    hdr["XRSLI13"] = 200.1
    hdr["YRSLI12"] = 300.1
    hdr["YRSLI13"] = 300.1

    hdr["XRSLI33"] = 1300.1
    hdr["YRSLI33"] = 1300.1
    hdr["XRSLI44"] = 1900.1
    hdr["YRSLI44"] = 1850.1
    hdr["XRSLI45"] = 1900.1
    hdr["YRSLI45"] = 1850.1
    #
    hdr["XVSLI12"] = 200.1
    hdr["XVSLI13"] = 200.1
    hdr["YVSLI12"] = 300.1
    hdr["YVSLI13"] = 300.1
    #
    hdr["XRSLI33"] = 1300.1
    hdr["YRSLI33"] = 1300.1
    hdr["XRSLI44"] = 1900.1
    hdr["YRSLI44"] = 1850.1
    hdr["XRSLI45"] = 1900.1
    hdr["YRSLI45"] = 1850.1
    return hdr


@pytest.mark.parametrize("hdr, nslits",[(create_test_header0(), 55), (create_test_header1(), 53)])
def test_csubar(hdr, nslits):
    dumdata = pkgutil.get_data('emirdrp.instrument.configs', 'bars_nominal_positions_test.txt')
    ss = S.StringIO(dumdata.decode('utf8'))
    bars_nominal_positions = numpy.loadtxt(ss)
    barmodel = create_bar_models(bars_nominal_positions)
    csu_conf = read_csu_from_header(barmodel, hdr)
    assert len(csu_conf.slits) == nslits


@pytest.mark.parametrize("hdr, nslits",[(create_test_header0(), 55), (create_test_header1(), 53)])
def test_merge_bars(hdr, nslits):
    mm = []
    for idx in range(1, EMIR_NBARS + 1):

        # References from header
        try:
            slit_t = hdr["SLIFL%d" % idx]
            target_type = TargetType(slit_t)
        except KeyError:
            target_type = TargetType.UNKNOWN

        xref = hdr.get("XRSLI%d" % idx, -100) - 1
        yref = hdr.get("YRSLI%d" % idx, -100) - 1
        target_coordinates = (xref, yref)

        xref = hdr.get("XVSLI%d" % idx, -100) - 1
        yref = hdr.get("YVSLI%d" % idx, -100) - 1
        target_coordinates_v = (xref, yref)

        mm.append((idx, target_type, target_coordinates, target_coordinates_v))

    bag = merge_slits(mm, max_slits=3, tol=1e-2)
    assert len(bag) == nslits


def test_csuconf1():
    dumdata = pkgutil.get_data('emirdrp.instrument.configs', 'bars_nominal_positions_test.txt')
    ss = S.StringIO(dumdata.decode('utf8'))
    bars_nominal_positions = numpy.loadtxt(ss)
    hdr = create_test_header1()
    barmodel = create_bar_models(bars_nominal_positions)
    csu_conf = read_csu_from_header(barmodel, hdr)
    assert csu_conf.is_open()
