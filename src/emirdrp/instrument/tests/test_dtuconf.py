
def compute_md5_checksum(conf):
    import hashlib
    mm = hashlib.md5()
    ndig = 3
    values = (conf.xdtu, conf.ydtu, conf.zdtu,
              conf.xdtu_0, conf.ydtu_0, conf.zdtu_0)
    m1 = [round(r, ndig) for r in values]
    m2 = ["{:+010.3f}".format(m) for m in m1]
    m3 = ':'.join(m2)
    mm.update(m3.encode('utf-8'))
    return mm.hexdigest()
