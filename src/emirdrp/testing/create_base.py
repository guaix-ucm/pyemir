import astropy.io.fits as fits
import numpy
from astropy.io import fits as fits


def compute_md5_checksum(conf):
    import hashlib

    mm = hashlib.md5()
    ndig = 3
    values = (conf.xdtu, conf.ydtu, conf.zdtu, conf.xdtu_0, conf.ydtu_0, conf.zdtu_0)
    m1 = [round(r, ndig) for r in values]
    m2 = ["{:+010.3f}".format(m) for m in m1]
    m3 = ":".join(m2)
    mm.update(m3.encode("utf-8"))
    return mm.hexdigest()


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


def create_frame_bpm():
    scene = numpy.zeros((10, 10), dtype="uint8")
    keys = {"UUID": "74263122-4747-4a68-b148-7834d702a733"}
    return create_image0(scene=scene, keys=keys)


def create_frame_dark():
    scene = numpy.zeros((10, 10))
    keys = {"UUID": "db4cf732-76a7-4f54-85bf-2db2f0ab46ef"}
    return create_image0(scene=scene, keys=keys)


def create_frame_flat():
    scene = numpy.zeros((10, 10))
    keys = {"UUID": "f7af05ee-27bd-4a04-957d-28d23b11a532"}
    return create_image0(scene=scene, keys=keys)


def create_images_mecs():
    aa = numpy.array([[1, 2, 3], [6, 5, 4]], dtype="float32")
    bb = numpy.array([[1]], dtype="uint16")
    uu = fits.PrimaryHDU(aa)
    mecs = fits.ImageHDU(bb, name="MECS")
    uu.header["TSUTC2"] = 0
    images = [fits.HDUList([uu, mecs]), fits.HDUList([uu, mecs])]
    return images
