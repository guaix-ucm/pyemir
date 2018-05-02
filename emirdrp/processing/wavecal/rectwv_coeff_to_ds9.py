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

from emirdrp.products import RectWaveCoeff

from emirdrp.core import EMIR_NAXIS1
from emirdrp.core import EMIR_NAXIS1_ENLARGED
from emirdrp.core import EMIR_NBARS

from numina.array.display.pause_debugplot import DEBUGPLOT_CODES


def save_ds9(output, filename):
    """Save ds9 region output info filename.

    Parameters
    ----------
    output : str
        String containing the full output to be exported as a ds9 region
        file.
    filename : str
        Output file name.

    """

    ds9_file = open(filename, 'wt')
    ds9_file.write(output)
    ds9_file.close()


def save_four_ds9(rectwv_coeff, debugplot=0):
    """Save the 4 possible ds9 region files.

    Parameters
    ----------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    debugplot : int
            Debugging level for messages and plots. For details see
            'numina.array.display.pause_debugplot.py'.

    """

    for limits, rectified, suffix in zip(
        ['frontiers', 'frontiers', 'boundaries', 'boundaries'],
        [False, True, False, True],
        ['rawimage', 'rectified', 'rawimage', 'rectified']
    ):
        output = rectwv_coeff_to_ds9(rectwv_coeff=rectwv_coeff,
                                     limits=limits,
                                     rectified=rectified)
        filename = 'ds9_' + limits + '_' + suffix + '.reg'
        if abs(debugplot) >= 10:
            print('>>> Saving: ', filename)
        save_ds9(output, filename)


def rectwv_coeff_to_ds9(rectwv_coeff,
                        limits=None,
                        rectified=False,
                        numpix=100):
    """Generate ds9 region output with requested slitlet limits

    Parameters
    ----------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    limits : str
        Region to be saved: the only two possibilities are 'boundaries'
        and 'frontiers'.
    rectified : bool
        If True, the regions correspond to the rectified image.
    numpix : int
        Number of points in which the X-range interval is subdivided
        in order to save each boundary as a connected set of line
        segments. This option is only relevant when rectified=False.

    Returns
    -------
    output : str
        String containing the full output to be exported as a ds9 region
        file.
    """

    # protections
    if limits not in ['boundaries', 'frontiers']:
        raise ValueError('Unexpect limits=' + str(limits))

    ds9_output = '# Region file format: DS9 version 4.1\n' \
                 'global color=green dashlist=2 4 width=2 ' \
                 'font="helvetica 10 normal roman" select=1 ' \
                 'highlite=1 dash=1 fixed=0 edit=1 ' \
                 'move=1 delete=1 include=1 source=1\nphysical\n#\n'

    ds9_output += '#\n# uuid (rectwv_coeff): {0}\n'.format(rectwv_coeff.uuid)
    ds9_output += \
        '# grism...............: {0}\n'.format(rectwv_coeff.tags['grism'])
    ds9_output += \
        '# filter..............: {0}\n'.format(rectwv_coeff.tags['filter'])

    for islitlet in range(1, EMIR_NBARS + 1):
        if islitlet not in rectwv_coeff.missing_slitlets:
            dumdict = rectwv_coeff.contents[islitlet - 1]
            if islitlet % 2 == 0:
                if limits == 'frontiers':
                    colorbox = '#0000ff'  # '#ff77ff'
                else:
                    colorbox = '#ff00ff'  # '#ff77ff'
            else:
                if limits == 'frontiers':
                    colorbox = '#0000ff'  # '#4444ff'
                else:
                    colorbox = '#00ffff'  # '#4444ff'

            ds9_output += '#\n# islitlet...........: {0}\n'.format(islitlet)
            ds9_output += '# csu_bar_slit_center: {0}\n'.format(
                dumdict['csu_bar_slit_center']
            )
            if rectified:
                if limits == 'frontiers':
                    ydum_lower = dumdict['y0_frontier_lower_expected']
                    ydum_upper = dumdict['y0_frontier_upper_expected']
                else:
                    ydum_lower = dumdict['y0_reference_lower_expected']
                    ydum_upper = dumdict['y0_reference_upper_expected']
                for ydum in [ydum_lower, ydum_upper]:
                    ds9_output += \
                        'line {0} {1} {2} {3}'.format(
                            1, ydum,
                            EMIR_NAXIS1_ENLARGED, ydum
                        )
                    ds9_output += ' # color={0}\n'.format(colorbox)
                ds9_output += 'text {0} {1} {{{2}}} # color={3} ' \
                              'font="helvetica 10 bold ' \
                              'roman"\n'.format(EMIR_NAXIS1_ENLARGED / 2 + 0.5,
                                                (ydum_lower + ydum_upper) / 2,
                                                islitlet,
                                                colorbox)
            else:
                if limits == 'frontiers':
                    pol_lower = Polynomial(
                        dumdict['frontier']['poly_coef_lower']
                    )
                    pol_upper = Polynomial(
                        dumdict['frontier']['poly_coef_upper']
                    )
                else:
                    pol_lower = Polynomial(
                        dumdict['spectrail']['poly_coef_lower']
                    )
                    pol_upper = Polynomial(
                        dumdict['spectrail']['poly_coef_upper']
                    )
                xdum = np.linspace(1, EMIR_NAXIS1, num=numpix)
                ydum = pol_lower(xdum)
                for i in range(len(xdum) - 1):
                    ds9_output += \
                        'line {0} {1} {2} {3}'.format(
                            xdum[i], ydum[i],
                            xdum[i + 1], ydum[i + 1]
                        )
                    ds9_output += ' # color={0}\n'.format(colorbox)
                ydum = pol_upper(xdum)
                for i in range(len(xdum) - 1):
                    ds9_output += \
                        'line {0} {1} {2} {3}'.format(
                            xdum[i], ydum[i],
                            xdum[i + 1], ydum[i + 1]
                        )
                    ds9_output += ' # color={0}\n'.format(colorbox)
                # slitlet label
                yc_lower = pol_lower(EMIR_NAXIS1 / 2 + 0.5)
                yc_upper = pol_upper(EMIR_NAXIS1 / 2 + 0.5)
                ds9_output += 'text {0} {1} {{{2}}} # color={3} ' \
                              'font="helvetica 10 bold ' \
                              'roman"\n'.format(EMIR_NAXIS1 / 2 + 0.5,
                                                (yc_lower + yc_upper) / 2,
                                                islitlet,
                                                colorbox)

    return ds9_output


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(
        description='description: generate ds9 region files associated '
                    'to a particular rectification and wavelength calibration'
    )

    # required arguments
    parser.add_argument("--rectwv_coeff", required=True,
                        help="Input JSON file with rectification and "
                             "wavelength calibration coefficients",
                        type=argparse.FileType('rt'))

    # optional arguments
    parser.add_argument("--debugplot",
                        help="Integer indicating plotting & debugging options"
                             " (default=0)",
                        default=0, type=int,
                        choices=DEBUGPLOT_CODES)
    parser.add_argument("--echo",
                        help="Display full command line",
                        action="store_true")
    args = parser.parse_args(args)

    if args.echo:
        print('\033[1m\033[31m% ' + ' '.join(sys.argv) + '\033[0m\n')

    # generate RectWaveCoeff object
    rectwv_coeff = RectWaveCoeff._datatype_load(
        args.rectwv_coeff.name)

    save_four_ds9(rectwv_coeff=rectwv_coeff, debugplot=args.debugplot)


if __name__ == "__main__":

    main()
