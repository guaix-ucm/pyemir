import astropy.io.fits as fits
import numpy


def create_image0(scene, hdr=None, keys=None):
    hdu = fits.PrimaryHDU(scene)
    if hdr:
        hdu.header = hdr
    if keys:
        for k, v in keys.items():
            hdu.header[k] = v

    return fits.HDUList([hdu])


def dither_pattern(center, base_angle, dist, npoints):

    base_angle_rad = base_angle / 180.0 * numpy.pi
    step = 2 * numpy.pi / npoints
    angles = base_angle_rad + numpy.arange(0.0, 2 * numpy.pi, step)
    x = center[0] + dist * numpy.cos(angles)
    y = center[1] + dist * numpy.sin(angles)
    return numpy.asarray([x, y]).T


def create_images_mecs():
    aa = numpy.array([[1, 2, 3], [6, 5, 4]], dtype="float32")
    bb = numpy.array([[1]], dtype="uint16")
    uu = fits.PrimaryHDU(aa)
    mecs = fits.ImageHDU(bb, name="MECS")
    uu.header["TSUTC2"] = 0
    images = [fits.HDUList([uu, mecs]), fits.HDUList([uu, mecs])]
    return images
