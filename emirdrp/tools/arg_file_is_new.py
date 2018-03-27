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

import os


def arg_file_is_new(parser, arg, mode='w'):
    """Auxiliary function to give an error if the file already exists.

    Parameters
    ----------
    parser : parser object
        Instance of argparse.ArgumentParser()
    arg : string
        File name.
    mode : string
        Optional string that specifies the mode in which the file is
        opened.

    Returns
    -------
    handler : file object
        Open file handle.

    """
    if os.path.exists(arg):
        parser.error("\nThe file\n\"%s\"\nalready exists and "
                     "cannot be overwritten!" % arg)
    else:
        # return an open file handle
        handler = open(arg, mode=mode)
        return handler
