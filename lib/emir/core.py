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

from numina.core.reciperesult import RecipeResultAutoQC as RecipeResult

def gather_info_dframe(dataframe):
    with dataframe.open() as hdulist:
        info = gather_info_hdu(hdulist)
    return info

def gather_info_hdu(hdulist):
    n_ext = len(hdulist)

    # READMODE is STRING
    readmode = hdulist[0].header.get('READMODE', 'undefined')
    bunit = hdulist[0].header.get('BUNIT', 'ADU')
    texp = hdulist[0].header.get('EXPTIME')
    adu_s = False

    if bunit.lower() == 'adu/s':
        adu_s = True
    

    return {'n_ext': n_ext, 
            'readmode': readmode,
            'texp': texp,
            'adu_s': adu_s}
    
def gather_info_frames(framelist):
    iinfo = []
    for frame in framelist:
        with frame.open() as hdulist:
            iinfo.append(gather_info_hdu(hdulist))
    return iinfo

EMIR_BIAS_MODES = ['simple']
