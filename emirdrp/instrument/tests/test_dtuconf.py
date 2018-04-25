
import pytest

from ..dtu_configuration import DtuConfiguration


@pytest.mark.xfail(reason="DtuConfiguration.__eq__ not complete")
def test_dtu_eq():

    dtu1 = DtuConfiguration()
    assert not (dtu1 == 100)


@pytest.mark.xfail(reason="DtuConfiguration.__ne__ not complete")
def test_dtu_eq():

    dtu1 = DtuConfiguration()
    assert dtu1 != 100


def compute_md5_checksum(conf):
    import hashlib
    mm = hashlib.md5()
    ndig = 3
    values = (conf.xdtu, conf.ydtu, conf.zdtu, conf.xdtu_0, conf.ydtu_0, conf.zdtu_0)
    m1 = [round(r, ndig) for r in values]
    m2 = ["{:+010.3f}".format(m) for m in m1]
    m3 = ':'.join(m2)
    mm.update(m3.encode('utf-8'))
    return mm.hexdigest()


def test_checksum():
    import itertools

    dtu1 = DtuConfiguration.define_from_values(
        -199.871994018555, -6.01818990707397, -535.546997070312,
        -200.02099609375, -5.93765020370483, -535.0
    )

    dtu2 = DtuConfiguration.define_from_values(
        -199.871994018555, -6.01818990707397, -535.546997070312,
        -200.02099609375, -5.93765020370483, -535.0
    )

    dtu3 = DtuConfiguration.define_from_values(
        -199.872, -6.01818990707397, -535.546997070312,
        -200.02099609375, -5.93765020370483, -535.0
    )

    dtu4 = DtuConfiguration.define_from_values(
        -100, -6, 535.546111,
        -1200.0, -5.937650, -535.0
    )

    for a, b in itertools.combinations([dtu1, dtu2, dtu3, dtu4], 2):
        ca = compute_md5_checksum(a)
        cb = compute_md5_checksum(b)
        if a == b:
            assert ca == cb
        else:
            assert ca != cb
