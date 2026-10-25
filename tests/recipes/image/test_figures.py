"""The figure of the image recipes is saved without pyplot"""

import types

import matplotlib.pyplot as plt
import numpy
from matplotlib.figure import Figure

from emirdrp.recipes.image.micro import MicroditheredImageRecipe


def test_recipe_figures(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    before = plt.get_fignums()
    recipe = MicroditheredImageRecipe(intermediate_results=True)
    assert isinstance(recipe._figure, Figure)
    data = numpy.random.default_rng(1).uniform(1.0, 2.0, size=(10, 10))
    recipe.figure_init(data.shape)
    recipe.figure_image(data, types.SimpleNamespace(lastname="frame"))
    recipe.figure_check_combination([1, 2, 3], [1.0, 1.1, 0.9], [3.0, 2.0, 1.5])
    recipe.figure_simple_image(data, title="image")
    recipe.figure_median_background([1.0, 1.2])
    recipe.figure_final_before_s(data)
    names = {path.name for path in tmp_path.iterdir()}
    assert names == {"figure-check-combination_i0.png", "figure-median-sky-background_i0.png"}
    assert plt.get_fignums() == before
