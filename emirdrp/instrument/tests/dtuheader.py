

# header fragment
def _dtu_header1():
    header = {}
    header['XDTU_F'] = 0.922918
    header['YDTU_F'] = 0.971815
    header['ZDTU_F'] = 1.0
    header['XDTU'] = -205.679000854492
    header['YDTU'] = -24.4878005981445
    header['ZDTU'] = -463.765991210938
    header['XDTU_0'] = -205.651000976562
    header['YDTU_0'] = -24.4794006347656
    header['ZDTU_0'] = -463.821014404297
    return header


def _dtu_header2():
    header = {}
    header['XDTU_F'] = 0.922917
    header['YDTU_F'] = 0.971815
    header['ZDTU_F'] = 1.0
    header['XDTU'] = -215.670
    header['YDTU'] = -241.487
    header['ZDTU'] = -261.765
    header['XDTU_0'] = -205.651
    header['YDTU_0'] = -21.4794
    header['ZDTU_0'] = -422.821
    return header


def _dtu_header3():
    import astropy.io.fits as fits
    header = fits.Header()
    header['XDTU_F'] = 0.922918
    header['YDTU_F'] = 0.971815
    header['ZDTU_F'] = 1.0
    header['XDTU'] = -225.670
    header['YDTU'] = -231.487
    header['ZDTU'] = -223.765
    header['XDTU_0'] = -295.651
    header['YDTU_0'] = -34.4794
    header['ZDTU_0'] = -432.821
    return header


HEADERS = [_dtu_header1(), _dtu_header2(), _dtu_header3()]