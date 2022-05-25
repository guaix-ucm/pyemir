
import pytest

from ..dtuconf import DtuConf, average, maxdiff
from ..dtu_configuration import DtuConfiguration, average_dtu_configurations
from ..dtu_configuration import maxdiff_dtu_configurations
from .dtuheader import HEADERS


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuc(hdr):
    dtu_a = DtuConf.from_header(hdr)
    dtu_b = DtuConfiguration.define_from_header(hdr)
    assert dtu_b.outdict() == dtu_a.outdict()


@pytest.mark.parametrize("hdr", HEADERS)
def test_dtuc_eq(hdr):
    dtu_a = DtuConf.from_header(hdr)
    dtu_b = DtuConfiguration.define_from_header(hdr)
    assert dtu_b == dtu_a
    assert dtu_a == dtu_b


@pytest.mark.parametrize("hdrs", [HEADERS])
def test_dtuc_neq(hdrs):
    import itertools
    for seq in itertools.combinations(HEADERS, 2):
        hdr1, hdr2 = seq
        dtu_a = DtuConf.from_header(hdr1)
        dtu_b = DtuConfiguration.define_from_header(hdr2)
        assert dtu_a != dtu_b
        assert dtu_b != dtu_a


@pytest.mark.parametrize("hdrs", [HEADERS])
def test_average(hdrs):
    confs1 = [DtuConf.from_header(hdr) for hdr in hdrs]
    confs2 = [DtuConfiguration.define_from_header(hdr) for hdr in hdrs]

    res1 = average(*confs1)
    res2 = average_dtu_configurations(confs2)

    assert res1.outdict() == res2.outdict()


@pytest.mark.parametrize("hdrs", [HEADERS])
def test_minmax(hdrs):
    confs1 = [DtuConf.from_header(hdr) for hdr in hdrs]
    confs2 = [DtuConfiguration.define_from_header(hdr) for hdr in hdrs]

    res1 = maxdiff(*confs1)
    res2 = maxdiff_dtu_configurations(confs2)

    assert res1.outdict() == res2.outdict()
