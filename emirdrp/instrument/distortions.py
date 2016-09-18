import numpy

# Platescale in radians
PLATESCALE_RADS = (0.1944 * 1 / 3600.0) * numpy. pi / 180.0


def exvp(pos_x, pos_y):
    pos_x = numpy.asarray(pos_x)
    pos_y = numpy.asarray(pos_y)
    # convert virtual pixel to real pixel
    # convert world coordinate to pixel
    center = [1024.5, 1024.5]
    cf  = PLATESCALE_RADS
    pos_base_x = pos_x - center[0]
    pos_base_y = pos_y - center[1]
    ra = numpy.hypot(pos_base_x, pos_base_y)
    thet = numpy.arctan2(pos_base_y, pos_base_x)

    r  = cf * ra

    rr1 = 1 + 14606.7 * r**2 + 1739716115.1 * r**4

    nx1 = rr1 * ra * numpy.cos(thet) + center[0]
    ny1 = rr1 * ra * numpy.sin(thet) + center[1]
    return nx1, ny1


def pvex(pos_x, pos_y):
    # convert real pixel to virtual pixel.
    # convert pixel to world coordinate
    pos_x = numpy.asarray(pos_x)
    pos_y = numpy.asarray(pos_y)

    center = [1024.5, 1024.5]
    cf  = PLATESCALE_RADS
    pos_base_x = pos_x - center[0]
    pos_base_y = pos_y - center[1]
    ra = numpy.hypot(pos_base_x, pos_base_y)
    thet = numpy.arctan2(pos_base_y, pos_base_x)

    r  = cf * ra
    rr1 = 1.000051 - 14892 * r**2 - 696254464 * r**4

    nx1 = rr1 * ra * numpy.cos(thet) + center[0]
    ny1 = rr1 * ra * numpy.sin(thet) + center[1]
    return nx1, ny1
