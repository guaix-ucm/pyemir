#
# Copyright 2008-2014 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyEmir is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
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
