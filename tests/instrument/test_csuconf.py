import io
import pkgutil

import numpy
import pytest

from emirdrp.instrument.csuconf import TargetType
from emirdrp.instrument.csuconf import create_bar_models, read_csu_from_header
from emirdrp.instrument.csuconf import merge_slits
from emirdrp.instrument.csuconf import EMIR_NBARS

from emirdrp.testing.create_headers import create_test_header0, create_test_header1


@pytest.mark.parametrize(
    "hdr, nslits", [(create_test_header0(), 55), (create_test_header1(), 53)]
)
def test_csubar(hdr, nslits):
    dumdata = pkgutil.get_data(
        "emirdrp.instrument.configs", "bars_nominal_positions_test.txt"
    )
    ss = io.StringIO(dumdata.decode("utf8"))
    bars_nominal_positions = numpy.loadtxt(ss)
    barmodel = create_bar_models(bars_nominal_positions)
    csu_conf = read_csu_from_header(barmodel, hdr)
    assert len(csu_conf.slits) == nslits


@pytest.mark.parametrize(
    "hdr, nslits", [(create_test_header0(), 55), (create_test_header1(), 53)]
)
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
    dumdata = pkgutil.get_data(
        "emirdrp.instrument.configs", "bars_nominal_positions_test.txt"
    )
    ss = io.StringIO(dumdata.decode("utf8"))
    bars_nominal_positions = numpy.loadtxt(ss)
    hdr = create_test_header1()
    barmodel = create_bar_models(bars_nominal_positions)
    csu_conf = read_csu_from_header(barmodel, hdr)
    assert csu_conf.is_open()
