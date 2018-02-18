from __future__ import division
from __future__ import print_function

import argparse
from copy import deepcopy
from datetime import datetime
import json
import os
import sys
from uuid import uuid4

from arg_file_is_new import arg_file_is_new

from emirdrp.core import EMIR_NBARS


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser()

    # positional arguments
    parser.add_argument("filename",
                        help="TXT with list of bounddict files",
                        type=argparse.FileType('r'))
    parser.add_argument("--outfile", required=True,
                        help="Output merged JSON file",
                        type=lambda x: arg_file_is_new(parser, x))

    # optional arguments
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")

    args = parser.parse_args()

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # initialize empty output
    bounddict = {}

    # read list of JSON files to be merged
    file_content = args.filename.read().splitlines()
    next_file_is_first = True
    for line in file_content:
        if len(line) > 0:
            if line[0] != '#':
                tmpfile = line.split()[0]
                if not os.path.isfile(tmpfile):
                    raise ValueError("File " + tmpfile + " not found!")
                tmpbounddict = json.loads(open(tmpfile).read())
                if next_file_is_first:
                    bounddict = deepcopy(tmpbounddict)
                    # update some values
                    bounddict['meta-info']['creation_date'] = \
                        datetime.now().isoformat()
                    bounddict['uuid'] = str(uuid4())
                    next_file_is_first = False
                else:
                    for islitlet in range(EMIR_NBARS):
                        cslitlet = "slitlet" + str(islitlet).zfill(2)
                        if cslitlet in tmpbounddict['contents']:
                            for dateobs in tmpbounddict['contents'][cslitlet]:
                                bounddict['contents'][cslitlet][dateobs] = \
                                    tmpbounddict['contents'][cslitlet][dateobs]

    # save merged JSON file
    with open(args.outfile.name, 'w') as fstream:
        json.dump(bounddict, fstream, indent=2, sort_keys=True)
        print('>>> Saving file ' + args.outfile.name)

if __name__ == "__main__":

    main()