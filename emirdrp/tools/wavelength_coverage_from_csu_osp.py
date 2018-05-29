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

import argparse
import numpy as np
from numpy.polynomial import Polynomial
import sys

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES
from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NBARS

def main(args=None):

    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: display arrangement of EMIR CSU bars'
    )

    # positional arguments
    parser.add_argument("csu_file",
                        help="TXT file with CSU configuration from the OSP",
                        type=argparse.FileType('rt'))
    parser.add_argument("configuration",
                        help="Instrument configuration (J, H, K, YJ, HK)",
                        choices=["J", "H", "K", "YJ", "HK"])

    # optional arguments
    parser.add_argument("--geometry",
                        help="Tuple x,y,dx,dy indicating window geometry",
                        default="0,0,640,480")
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting & debugging options"
                             " (default=12)",
                        default=12, type=int,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")
    args = parser.parse_args(args)

    if args.echo:
        print('\033[1m\033[31mExecuting: ' + ' '.join(sys.argv) + '\033[0m\n')

    # geometry
    if args.geometry is None:
        geometry = None
    else:
        tmp_str = args.geometry.split(",")
        x_geom = int(tmp_str[0])
        y_geom = int(tmp_str[1])
        dx_geom = int(tmp_str[2])
        dy_geom = int(tmp_str[3])
        geometry = x_geom, y_geom, dx_geom, dy_geom

    # initialize empty array to store location of bars
    csu_bar_left = np.zeros(EMIR_NBARS)
    csu_bar_right = np.zeros(EMIR_NBARS)
    csu_bar_slit_center = np.zeros(EMIR_NBARS)
    csu_bar_slit_width = np.zeros(EMIR_NBARS)

    # initialize empty array to store flag
    csu_flag = np.zeros(2*EMIR_NBARS, dtype=int)

    # determine calibration
    if args.configuration == "J":
        crval1 = Polynomial([1.25137094e+04, -4.81553731e+00, 4.70039758e-04])
        cdelt1 = Polynomial([7.74133267e-01, -4.72423718e-05, 2.79842624e-08])
    elif args.configuration == "H":
        crval1 = Polynomial([1.65536274e+04, -7.63517173e+00, 7.74790265e-04])
        cdelt1 = Polynomial([1.21327515e+00, 1.42140078e-05, -1.27489119e-07])
    elif args.configuration == "K":
        crval1 = Polynomial([2.21044741e+04, -1.08737529e+01, 9.05081653e-04])
        cdelt1 = Polynomial([1.72696857e+00, 2.35009351e-05, -1.02164228e-07])
    elif args.configuration == "YJ":
        crval1 = Polynomial([1.04272465e+04, -2.33176855e+01, 6.55101267e-03])
        cdelt1 = Polynomial([3.49037727e+00, 1.26008332e-03, -4.66149678e-06])
    elif args.configuration == "HK":
        crval1 = Polynomial([2.00704978e+04, -4.07702886e+01, -5.95247468e-03])
        cdelt1 = Polynomial([6.54247758e+00, 2.09061196e-03, -2.48206609e-06])
    else:
        raise ValueError("Invalid configuration: " + args.configuration)

    # read input file
    file_content = args.csu_file.read().splitlines()
    next_id_bar = 1
    for line in file_content:
        if len(line) > 0:
            if line[0] not in ['#']:
                line_contents = line.split()
                id_bar = int(line_contents[0])
                position = float(line_contents[1])
                flag = int(line_contents[2])
                if id_bar == next_id_bar:
                    csu_flag[id_bar-1] = flag
                    if id_bar <= EMIR_NBARS:
                        csu_bar_left[id_bar-1] = position
                        next_id_bar = id_bar + EMIR_NBARS
                    else:
                        csu_bar_right[id_bar - EMIR_NBARS -1] = \
                                341.5 - position
                        next_id_bar = id_bar - EMIR_NBARS + 1
                else:
                    raise ValueError("Unexpected id_bar:" + str(id_bar))

    # compute csu_bar_slit_center
    print("bar    left   right  center   width flag   min.wave  max.wave")
    print("--- ------- ------- -------   ----- ----   --------  --------")
    for i in range(EMIR_NBARS):
        csu_bar_slit_center[i] = (csu_bar_left[i] + csu_bar_right[i])/2
        csu_bar_slit_width[i] = csu_bar_right[i] - csu_bar_left[i]
        csu_crval1 = crval1(csu_bar_slit_center[i])
        csu_cdelt1 = cdelt1(csu_bar_slit_center[i])
        csu_crvaln = csu_crval1 + (EMIR_NAXIS1 - 1) * csu_cdelt1
        print(" {0:2d} {1:7.3f} {2:7.3f} {3:7.3f} {4:7.3f}   {5:2d}  "
              " {6:8.2f}  {7:8.2f}".format(
              i+1,
              csu_bar_left[i], csu_bar_right[i],
              csu_bar_slit_center[i], csu_bar_slit_width[i], csu_flag[i],
              csu_crval1, csu_crvaln)
        )


if __name__ == "__main__":
    main()