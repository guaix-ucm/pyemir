from __future__ import division
from __future__ import print_function

import os


def arg_file_is_new(parser, arg):
    """Auxiliary function to give an error if the file already exists.

    Parameters
    ----------
    parser : parser object
        Instance of argparse.ArgumentParser()
    arg : string
        File name.

    Returns
    -------
    handler : file object
        Open file handle.

    """
    if os.path.exists(arg):
        parser.error("\nThe file\n\"%s\"\nalready exist and "
                     "cannot be overwritten!" % arg)
    else:
        # return an open file handle
        handler = open(arg, 'w')
        return handler
