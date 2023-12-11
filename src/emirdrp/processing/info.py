#
# Copyright 2008-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#


def gather_info_dframe(dataframe):
    with dataframe.open() as hdulist:
        info = gather_info_hdu(hdulist)
    return info


def gather_info_hdu(hdulist):

    _meta = {'readmode': ('READMODE', 'undefined'),
             'texp': ('EXPTIME', None),
             'grism': ('GRISM', 'undefined'),
             'filter': ('FILTER', 'undefined'),
             'obsmode': ('OBSMODE', 'undefined'),
             'tstamp': ('TSTAMP', 'undefined')
             }

    # READMODE is STRING
    meta = {}
    meta['n_ext'] = len(hdulist)
    extnames = [hdu.header.get('extname', '') for hdu in hdulist[1:]]
    meta['name_ext'] = ['PRIMARY'] + extnames
    for key, val in _meta.items():
        meta[key] = hdulist[0].header.get(val[0], val[1])

    return meta


def gather_info_frames(framelist):
    iinfo = []
    for frame in framelist:
        with frame.open() as hdulist:
            iinfo.append(gather_info_hdu(hdulist))
    return iinfo


def gather_info(recipeinput):
    from numina.core import DataFrame
    klass = recipeinput.__class__
    metadata = {}
    for key in klass.stored():
        val = getattr(recipeinput, key)
        if isinstance(val, DataFrame):
            metadata[key] = gather_info_dframe(val)
        elif hasattr(val, 'frames'):
            metas = []
            for f in val.frames:
                metas.append(gather_info_dframe(f))
            metadata[key] = metas
        else:
            pass
    return metadata
