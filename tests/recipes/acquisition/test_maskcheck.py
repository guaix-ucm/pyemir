"""The debug plot of the centroid is drawn without pyplot"""

import matplotlib.pyplot as plt
import numpy
from matplotlib.figure import Figure
from numina.array.bbox import BoundingBox

from emirdrp.recipes.acquisition.maskcheck import comp_centroid


def test_comp_centroid_debug_plot():
    y, x = numpy.indices((100, 100))
    data = numpy.random.default_rng(1).normal(100.0, 1.0, size=(100, 100))
    data += 1000 * numpy.exp(-0.5 * ((x - 52) ** 2 + (y - 47) ** 2) / 2.0**2)
    before = plt.get_fignums()
    xc, yc, ax = comp_centroid(data, BoundingBox(30, 70, 30, 70), debug_plot=True, plot_reference=(50, 50))
    # the centroid is relative to the bounding box
    assert abs(xc + 30 - 52) < 0.1 and abs(yc + 30 - 47) < 0.1
    assert isinstance(ax.figure, Figure)
    assert plt.get_fignums() == before
