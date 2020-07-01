from numina.array import combine
import numpy

from ..combine import combine_images, scale_with_median


def create_images():
    import astropy.io.fits as fits
    aa = numpy.array([[1, 2, 3], [6, 5, 4]], dtype='float32')
    bb = numpy.array([[1]], dtype='uint16')
    uu = fits.PrimaryHDU(aa)
    mecs = fits.ImageHDU(bb, name='MECS')
    uu.header['TSUTC2'] = 0
    images = [fits.HDUList([uu, mecs]), fits.HDUList([uu, mecs])]
    return images


def test_combine_images():
    images = create_images()
    res = combine_images(images)
    assert 'MECS' in res


def test_combine_images2():
    images = create_images()
    res = combine_images(images, errors=True)
    assert 'MECS' in res
    assert 'VARIANCE' in res
    assert 'MAP' in res


def test_combine_images3():
    images = create_images()
    scaled_mean = scale_with_median(combine.mean)
    res = combine_images(images, method=scaled_mean, errors=True)
    assert 'MECS' in res
    assert 'VARIANCE' in res
    assert 'MAP' in res
