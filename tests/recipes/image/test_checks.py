"""The figure of the photometry check is saved without pyplot"""

import types

import matplotlib.pyplot as plt
import numpy
from astropy.io import fits

from emirdrp.recipes.image.checks import check_photometry
from emirdrp.recipes.image.naming import name_skysub_proc

SHAPE = (200, 200)
STARS = [(80, 90), (110, 100), (95, 125), (120, 80)]


def star_image(rng):
    y, x = numpy.indices(SHAPE)
    data = rng.normal(100.0, 1.0, size=SHAPE)
    for xc, yc in STARS:
        data += 1000 * numpy.exp(-0.5 * ((x - xc) ** 2 + (y - yc) ** 2) / 2.0**2)
    return data.astype("float32")


def test_check_photometry_figure(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rng = numpy.random.default_rng(1)
    fits.writeto("result_i0.fits", star_image(rng))
    frames = []
    for idx in range(3):
        frame = types.SimpleNamespace(baselabel=f"frame{idx}", airmass=1.0, valid_science=True)
        header = fits.Header({"AIRMASS": 1.0, "TSTAMP": 10.0 * idx})
        fits.writeto(name_skysub_proc(frame.baselabel, 0), star_image(rng), header)
        frames.append(frame)
    before = plt.get_fignums()
    check_photometry(frames, [numpy.zeros(SHAPE)], seeing_fwhm=None, border=50, save_figure=True)
    assert (tmp_path / "figure-relative-flux_i0.png").read_bytes()[:4] == b"\x89PNG"
    assert plt.get_fignums() == before
