from numina.array import combine

from emirdrp.processing.combine import combine_images, scale_with_median
from emirdrp.testing.create_base import create_images_mecs


def test_combine_images():
    images = create_images_mecs()
    res = combine_images(images)
    assert "MECS" in res


def test_combine_images2():
    images = create_images_mecs()
    res = combine_images(images, errors=True)
    assert "MECS" in res
    assert "VARIANCE" in res
    assert "MAP" in res


def test_combine_images3():
    images = create_images_mecs()
    scaled_mean = scale_with_median(combine.mean)
    res = combine_images(images, method=scaled_mean, errors=True)
    assert "MECS" in res
    assert "VARIANCE" in res
    assert "MAP" in res
