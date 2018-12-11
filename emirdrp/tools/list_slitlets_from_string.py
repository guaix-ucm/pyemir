#
# Copyright 2008-2018 Universidad Complutense de Madrid
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

from __future__ import division
from __future__ import print_function

import re


def list_slitlets_from_string(s, islitlet_min, islitlet_max):
    """Return list of slitlets from string specification.

    Parameters
    ----------
    s : string
        String defining the slitlets. The slitlets must be specify
        as a set of n1[,n2[,step]] tuples. If only n1 is provided,
        the slitlet number n1 is considered. When n1 and n2 are given
        but step is missing, step=1 is assumed. Finally, when all
        n1, n2 and step are given, slitlets considered are those
        returned by range(n1, n2 + 1, step) function.
    islitlet_min : int
        Minimum slitlet number allowed.
    islitlet_max : int
        Maximum slitlet number allowed.

    Returns
    -------
    list_slitlets : Python list
        List of slitlets.

    """

    # protection
    if not isinstance(s, str):
        print('type(s): ', type(s))
        print('ERROR: function expected a string parameter')

    # initialize empty output
    set_slitlets = set()

    # remove leading blank spaces
    s = re.sub('^ *', '', s)
    # remove trailing blank spaces
    s = re.sub(' *$', '', s)
    # remove blank space before ','
    s = re.sub(' *,', ',', s)
    # remove blank spaces after  ','
    s = re.sub(', *', ',', s)
    # remove many blank spaces by a single blank space
    s = re.sub(' +', ' ', s)
    stuples = s.split()
    for item in stuples:
        subitems = item.split(',')
        nsubitems = len(subitems)
        if nsubitems == 1:
            n1 = int(subitems[0])
            n2 = n1
            step = 1
        elif nsubitems == 2:
            n1 = int(subitems[0])
            n2 = int(subitems[1])
            step = 1
        elif nsubitems == 3:
            n1 = int(subitems[0])
            n2 = int(subitems[1])
            step = int(subitems[2])
        else:
            raise ValueError('Unexpected slitlet range:', s)
        for i in range(n1, n2 + 1, step):
            if islitlet_min <= i <= islitlet_max:
                set_slitlets.add(i)
            else:
                print('islitlet_min: ', islitlet_min)
                print('islitlet_max: ', islitlet_max)
                print('i...........: ', i)
                raise ValueError("Slitlet number out of range!")

    list_slitlets = list(set_slitlets)
    list_slitlets.sort()

    # return result
    return list_slitlets
