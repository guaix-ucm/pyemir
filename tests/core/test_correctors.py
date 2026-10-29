"""The getters of the correctors of the calibrations"""

import types

import astropy.io.fits as fits
import numina.processing as proc
import numpy
import pytest
from numina.core import DataFrame
from numina.util.node import IdNode

import emirdrp.core.correctors as cor
from emirdrp.datamodel import EmirDataModel
from emirdrp.processing.flatfield import FlatFieldCorrector

UUID = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"


def calibration(value=2.0):
    hdu = fits.PrimaryHDU(numpy.full((10, 20), value, dtype="float32"))
    hdu.header["UUID"] = UUID
    return DataFrame(frame=fits.HDUList([hdu]))


# key of the requirement, getter, class of the corrector
CORRECTORS = [
    ("master_bpm", cor.get_corrector_p, proc.BadPixelCorrector),
    ("master_sky", cor.get_corrector_s, proc.SkyCorrector),
    ("master_dark", cor.get_corrector_d, proc.DarkCorrector),
    ("master_flat", cor.get_corrector_f, FlatFieldCorrector),
    ("master_flat", cor.get_corrector_sf, FlatFieldCorrector),
]


@pytest.mark.parametrize("key, getter, corrector_class", CORRECTORS)
def test_get_corrector(key, getter, corrector_class):
    """The corrector of the calibration, with its id and the kind of id"""
    rinput = types.SimpleNamespace(**{key: calibration()})
    corrector = getter(rinput, {key: {}}, None, EmirDataModel())
    assert isinstance(corrector, corrector_class)
    assert corrector.calibid == f"uuid:{UUID}"


@pytest.mark.parametrize("key, getter, corrector_class", CORRECTORS)
def test_get_corrector_missing(key, getter, corrector_class):
    """Without the calibration, the corrector does nothing"""
    rinput = types.SimpleNamespace(**{key: None})
    assert isinstance(getter(rinput, {}, None, EmirDataModel()), IdNode)


@pytest.mark.parametrize("readmode, corrector_class", [("simple", proc.BiasCorrector), ("cds", IdNode)])
def test_get_corrector_bias(readmode, corrector_class):
    """The bias is only used with some read modes"""
    rinput = types.SimpleNamespace(master_bias=calibration())
    meta = {"obresult": [{"readmode": readmode}]}
    corrector = cor.get_corrector_b(rinput, meta, None, EmirDataModel())
    assert isinstance(corrector, corrector_class)


def test_get_corrector_bias_without_images():
    rinput = types.SimpleNamespace(master_bias=calibration())
    assert isinstance(cor.get_corrector_b(rinput, {"obresult": []}, None, EmirDataModel()), IdNode)


def test_get_corrector_gen_deprecated():
    rinput = types.SimpleNamespace(master_dark=calibration())
    with pytest.warns(DeprecationWarning, match="use numina.processing.getters.load_corrector"):
        corrector = cor.get_corrector_gen(rinput, EmirDataModel(), proc.DarkCorrector, "master_dark")
    assert isinstance(corrector, proc.DarkCorrector)
