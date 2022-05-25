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
import pkgutil
from six import StringIO
import sys

from emirdrp.processing.wavecal.set_wv_parameters import set_wv_parameters
from emirdrp.products import RectWaveCoeff

from numina.array.distortion import fmap

from emirdrp.core import EMIR_NAXIS1
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

    # retrieve relevant wavelength calibration parameters
    grism_name = rectwv_coeff.tags['grism']
    filter_name = rectwv_coeff.tags['filter']
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    naxis1_enlarged = wv_parameters['naxis1_enlarged']
    crpix1_enlarged = wv_parameters['crpix1_enlarged']
    crval1_enlarged = wv_parameters['crval1_enlarged']
    cdelt1_enlarged = wv_parameters['cdelt1_enlarged']

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
                crpix1_linear = 1.0
                crval1_linear = dumdict['crval1_linear']
                cdelt1_linear = dumdict['cdelt1_linear']
                if limits == 'frontiers':
                    ydum_lower = dumdict['y0_frontier_lower_expected']
                    ydum_upper = dumdict['y0_frontier_upper_expected']
                else:
                    ydum_lower = dumdict['y0_reference_lower_expected']
                    ydum_upper = dumdict['y0_reference_upper_expected']
                wave_ini = crval1_linear + \
                           (0.5 - crpix1_linear) * cdelt1_linear
                xdum_ini = (wave_ini - crval1_enlarged) / cdelt1_enlarged
                xdum_ini += crpix1_enlarged
                wave_end = crval1_linear + \
                           (EMIR_NAXIS1 + 0.5 - crpix1_linear) * cdelt1_linear
                xdum_end = (wave_end - crval1_enlarged) / cdelt1_enlarged
                xdum_end += crpix1_enlarged
                for ydum in [ydum_lower, ydum_upper]:
                    ds9_output += \
                        'line {0} {1} {2} {3}'.format(
                            xdum_ini, ydum,
                            xdum_end, ydum
                        )
                    ds9_output += ' # color={0}\n'.format(colorbox)
                # slitlet label
                ydum_label = (ydum_lower + ydum_upper) / 2.0
                xdum_label = EMIR_NAXIS1 / 2 + 0.5
                wave_center = crval1_linear + \
                              (xdum_label - crpix1_linear) * cdelt1_linear
                xdum_label = (wave_center - crval1_enlarged) / cdelt1_enlarged
                xdum_label += crpix1_enlarged
                ds9_output += 'text {0} {1} {{{2}}} # color={3} ' \
                              'font="helvetica 10 bold ' \
                              'roman"\n'.format(xdum_label, ydum_label,
                                                islitlet, colorbox)
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
                ydum -= float(rectwv_coeff.global_integer_offset_y_pix)
                for i in range(len(xdum) - 1):
                    ds9_output += \
                        'line {0} {1} {2} {3}'.format(
                            xdum[i], ydum[i],
                            xdum[i + 1], ydum[i + 1]
                        )
                    ds9_output += ' # color={0}\n'.format(colorbox)
                ydum = pol_upper(xdum)
                ydum -= float(rectwv_coeff.global_integer_offset_y_pix)
                for i in range(len(xdum) - 1):
                    ds9_output += \
                        'line {0} {1} {2} {3}'.format(
                            xdum[i], ydum[i],
                            xdum[i + 1], ydum[i + 1]
                        )
                    ds9_output += ' # color={0}\n'.format(colorbox)
                # slitlet label
                xdum_label = EMIR_NAXIS1 / 2 + 0.5
                ydum_lower = pol_lower(xdum_label)
                ydum_upper = pol_upper(xdum_label)
                ydum_label = (ydum_lower + ydum_upper) / 2.0
                ds9_output += 'text {0} {1} {{{2}}} # color={3} ' \
                              'font="helvetica 10 bold ' \
                              'roman"\n'.format(xdum_label, ydum_label,
                                                islitlet, colorbox)

    return ds9_output


def save_spectral_lines_ds9(rectwv_coeff, debugplot=0):
    """Save expected location of arc and OH airglow to ds9 region files.

    Parameters
    ----------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    debugplot : int
            Debugging level for messages and plots. For details see
            'numina.array.display.pause_debugplot.py'.

    """

    for spectral_lines, rectified, suffix in zip(
        ['arc', 'arc', 'oh', 'oh'],
        [False, True, False, True],
        ['rawimage', 'rectified', 'rawimage', 'rectified']
    ):
        output = spectral_lines_to_ds9(rectwv_coeff=rectwv_coeff,
                                       spectral_lines=spectral_lines,
                                       rectified=rectified)
        filename = 'ds9_' + spectral_lines + '_' + suffix + '.reg'
        if abs(debugplot) >= 10:
            print('>>> Saving: ', filename)
        save_ds9(output, filename)


def spectral_lines_to_ds9(rectwv_coeff,
                    spectral_lines=None,
                    rectified=False):
    """Generate ds9 region output with requested spectral lines

    Parameters
    ----------
    rectwv_coeff : RectWaveCoeff instance
        Rectification and wavelength calibration coefficients for the
        particular CSU configuration.
    spectral_lines : str
        Spectral lines to be saved: the only two possibilities are 'arc'
        and 'oh'.
    rectified : bool
        If True, the regions correspond to the rectified image.

    Returns
    -------
    output : str
        String containing the full output to be exported as a ds9 region
        file.
    """

    # protections
    if spectral_lines not in ['arc', 'oh']:
        raise ValueError('Unexpected spectral lines=' + str(spectral_lines))

    if spectral_lines == 'arc':
        grism_name = rectwv_coeff.tags['grism']
        if grism_name == 'LR':
            catlines_file = 'lines_argon_neon_xenon_empirical_LR.dat'
        else:
            catlines_file = 'lines_argon_neon_xenon_empirical.dat'
        dumdata = pkgutil.get_data('emirdrp.instrument.configs',
                                   catlines_file)
        arc_lines_tmpfile = StringIO(dumdata.decode('utf8'))
        catlines = np.genfromtxt(arc_lines_tmpfile)
        # define wavelength and flux as separate arrays
        catlines_all_wave = catlines[:, 0]
        catlines_all_flux = catlines[:, 1]
    elif spectral_lines == 'oh':
        dumdata = pkgutil.get_data('emirdrp.instrument.configs',
                                   'Oliva_etal_2013.dat')
        oh_lines_tmpfile = StringIO(dumdata.decode('utf8'))
        catlines = np.genfromtxt(oh_lines_tmpfile)
        # define wavelength and flux as separate arrays
        catlines_all_wave = np.concatenate((catlines[:, 1],
                                            catlines[:, 0]))
        catlines_all_flux = np.concatenate((catlines[:, 2],
                                            catlines[:, 2]))
    else:
        raise ValueError('This should not happen!')

    # retrieve relevant wavelength calibration parameters
    grism_name = rectwv_coeff.tags['grism']
    filter_name = rectwv_coeff.tags['filter']
    wv_parameters = set_wv_parameters(filter_name, grism_name)
    naxis1_enlarged = wv_parameters['naxis1_enlarged']
    crpix1_enlarged = wv_parameters['crpix1_enlarged']
    crval1_enlarged = wv_parameters['crval1_enlarged']
    cdelt1_enlarged = wv_parameters['cdelt1_enlarged']

    ds9_output = '# Region file format: DS9 version 4.1\n' \
                 'global color=#00ffff dashlist=0 0 width=2 ' \
                 'font="helvetica 10 normal roman" select=1 ' \
                 'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 ' \
                 'include=1 source=1\nphysical\n'

    ds9_output += '#\n# uuid (rectwv_coeff): {0}\n'.format(
        rectwv_coeff.uuid)
    ds9_output += \
        '# grism...............: {0}\n'.format(rectwv_coeff.tags['grism'])
    ds9_output += \
        '# filter..............: {0}\n'.format(rectwv_coeff.tags['filter'])

    global_integer_offset_x_pix = rectwv_coeff.global_integer_offset_x_pix
    global_integer_offset_y_pix = rectwv_coeff.global_integer_offset_y_pix

    for islitlet in range(1, EMIR_NBARS + 1):
        if islitlet not in rectwv_coeff.missing_slitlets:
            dumdict = rectwv_coeff.contents[islitlet - 1]
            if islitlet % 2 == 0:
                colorbox = '#ff00ff'  # '#ff77ff'
            else:
                colorbox = '#00ffff'  # '#4444ff'

            ds9_output += '#\n# islitlet...........: {0}\n'.format(
                islitlet)
            ds9_output += '# csu_bar_slit_center: {0}\n'.format(
                dumdict['csu_bar_slit_center']
            )
            crpix1_linear = 1.0
            crval1_linear = dumdict['crval1_linear']
            cdelt1_linear = dumdict['cdelt1_linear']
            wave_ini = crval1_linear + \
                       (0.5 - crpix1_linear) * cdelt1_linear
            wave_end = crval1_linear + \
                       (EMIR_NAXIS1 + 0.5 - crpix1_linear) * cdelt1_linear
            if rectified:
                ydum_lower = dumdict['y0_reference_lower_expected']
                ydum_upper = dumdict['y0_reference_upper_expected']
                # spectral lines
                for wave in catlines_all_wave:
                    if wave_ini <= wave <= wave_end:
                        xdum = (wave - crval1_enlarged) / cdelt1_enlarged
                        xdum += crpix1_enlarged
                        ds9_output += \
                            'line {0} {1} {2} {3}'.format(
                                xdum, ydum_lower,
                                xdum, ydum_upper
                            )
                        ds9_output += ' # color={0}\n'.format(colorbox)
                # slitlet label
                ydum_label = (ydum_lower + ydum_upper) / 2.0
                xdum_label = EMIR_NAXIS1 / 2 + 0.5
                wave_center = crval1_linear + \
                              (xdum_label - crpix1_linear) * cdelt1_linear
                xdum_label = (wave_center - crval1_enlarged) / cdelt1_enlarged
                xdum_label += crpix1_enlarged
                ds9_output += 'text {0} {1} {{{2}}} # color={3} ' \
                              'font="helvetica 10 bold ' \
                              'roman"\n'.format(xdum_label, ydum_label,
                                                islitlet, colorbox)
            else:
                bb_ns1_orig = dumdict['bb_ns1_orig']
                ttd_order = dumdict['ttd_order']
                aij = dumdict['ttd_aij']
                bij = dumdict['ttd_bij']
                min_row_rectified = float(dumdict['min_row_rectified'])
                max_row_rectified = float(dumdict['max_row_rectified'])
                mean_row_rectified = (min_row_rectified + max_row_rectified)/2
                wpoly_coeff = dumdict['wpoly_coeff']
                x0 = []
                y0 = []
                x1 = []
                y1 = []
                x2 = []
                y2 = []
                # spectral lines
                for wave in catlines_all_wave:
                    if wave_ini <= wave <= wave_end:
                        tmp_coeff = np.copy(wpoly_coeff)
                        tmp_coeff[0] -= wave
                        tmp_xroots = np.polynomial.Polynomial(
                            tmp_coeff).roots()
                        for dum in tmp_xroots:
                            if np.isreal(dum):
                                dum = dum.real
                                if 1 <= dum <= EMIR_NAXIS1:
                                    x0.append(dum)
                                    y0.append(mean_row_rectified)
                                    x1.append(dum)
                                    y1.append(min_row_rectified)
                                    x2.append(dum)
                                    y2.append(max_row_rectified)
                        pass
                xx0, yy0 = fmap(ttd_order, aij, bij, np.array(x0),
                                np.array(y0))
                xx0 -= global_integer_offset_x_pix
                yy0 += bb_ns1_orig
                yy0 -= global_integer_offset_y_pix
                xx1, yy1 = fmap(ttd_order, aij, bij, np.array(x1),
                                np.array(y1))
                xx1 -= global_integer_offset_x_pix
                yy1 += bb_ns1_orig
                yy1 -= global_integer_offset_y_pix
                xx2, yy2 = fmap(ttd_order, aij, bij, np.array(x2),
                                np.array(y2))
                xx2 -= global_integer_offset_x_pix
                yy2 += bb_ns1_orig
                yy2 -= global_integer_offset_y_pix
                for xx1_, xx2_, yy1_, yy2_ in zip(xx1, xx2, yy1, yy2):
                    ds9_output += \
                        'line {0} {1} {2} {3}'.format(
                            xx1_, yy1_, xx2_, yy2_
                        )
                    ds9_output += ' # color={0}\n'.format(colorbox)
                # slitlet label
                pol_lower = Polynomial(dumdict['spectrail']['poly_coef_lower'])
                pol_upper = Polynomial(dumdict['spectrail']['poly_coef_upper'])
                xdum_label = EMIR_NAXIS1 / 2 + 0.5
                ydum_lower = pol_lower(xdum_label)
                ydum_upper = pol_upper(xdum_label)
                ydum_label = (ydum_lower + ydum_upper) / 2.0
                ds9_output += 'text {0} {1} {{{2}}} # color={3} ' \
                              'font="helvetica 10 bold ' \
                              'roman"\n'.format(xdum_label, ydum_label,
                                                islitlet, colorbox)

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
    save_spectral_lines_ds9(rectwv_coeff=rectwv_coeff,
                            debugplot=args.debugplot)


if __name__ == "__main__":

    main()
