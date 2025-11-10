def create_test_header0():
    hdr = {}
    for i in range(55):
        hdr["CSUP{}".format(i + 1)] = -100
    for i in range(55, 110):
        hdr["CSUP{}".format(i + 1)] = -100
    return hdr


def create_test_header1():
    hdr = create_test_header0()

    hdr["SLIFL12"] = 2
    hdr["SLIFL13"] = 2
    hdr["SLIFL33"] = 2
    hdr["SLIFL44"] = 2
    hdr["SLIFL45"] = 2

    hdr["XRSLI12"] = 200.1
    hdr["XRSLI13"] = 200.1
    hdr["YRSLI12"] = 300.1
    hdr["YRSLI13"] = 300.1

    hdr["XRSLI33"] = 1300.1
    hdr["YRSLI33"] = 1300.1
    hdr["XRSLI44"] = 1900.1
    hdr["YRSLI44"] = 1850.1
    hdr["XRSLI45"] = 1900.1
    hdr["YRSLI45"] = 1850.1
    #
    hdr["XVSLI12"] = 200.1
    hdr["XVSLI13"] = 200.1
    hdr["YVSLI12"] = 300.1
    hdr["YVSLI13"] = 300.1
    #
    hdr["XRSLI33"] = 1300.1
    hdr["YRSLI33"] = 1300.1
    hdr["XRSLI44"] = 1900.1
    hdr["YRSLI44"] = 1850.1
    hdr["XRSLI45"] = 1900.1
    hdr["YRSLI45"] = 1850.1
    return hdr


def create_dtu_header_example():
    header = dict()
    header["XDTU_F"] = 0.922918
    header["YDTU_F"] = 0.971815
    header["ZDTU_F"] = 1.0
    header["XDTU"] = -205.679000854492
    header["YDTU"] = -24.4878005981445
    header["ZDTU"] = -463.765991210938
    header["XDTU_0"] = -205.651000976562
    header["YDTU_0"] = -24.4794006347656
    header["ZDTU_0"] = -463.821014404297
    return header
