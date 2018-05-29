#
# Copyright 2011-2018 Universidad Complutense de Madrid
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

"""Routines shared by image mode recipes"""

import logging
import operator

import six
import numpy
from astropy.io import fits
import sep
import matplotlib.pyplot as plt

from emirdrp.util.sextractor import SExtractor
from .naming import name_skysub_proc


_logger = logging.getLogger(__name__)

# Actions to carry over images when checking the flux
# of the objects in different images


def warn_action(img):
    _logger.warning('Image %s has low flux in objects', img.baselabel)
    img.valid_science = True


def reject_action(img):
    img.valid_science = False
    _logger.info('Image %s rejected, has low flux in objects', img.baselabel)
    pass


def default_action(img):
    _logger.info(
        'Image %s accepted, has correct flux in objects', img.baselabel)
    img.valid_science = True

# Actions
_dactions = {'warn': warn_action,
             'reject': reject_action, 'default': default_action}


def check_photometry(frames, sf_data, seeing_fwhm, step=0,
                     border=300, extinction=0.0,
                     check_photometry_levels=[0.5, 0.8],
                     check_photometry_actions=['warn', 'warn', 'default'],
                     figure=None):
    # Check photometry of few objects
    weigthmap = 'weights4rms.fits'

    wmap = numpy.ones_like(sf_data[0], dtype='bool')

    # Center of the image
    wmap[border:-border, border:-border] = 0
    # fits.writeto(weigthmap, wmap.astype('uintt8'), overwrite=True)

    basename = 'result_i%0d.fits' % (step)

    data_res = fits.getdata(basename)
    data_res = data_res.byteswap().newbyteorder()
    bkg = sep.Background(data_res)
    data_sub = data_res - bkg

    _logger.info('Runing source extraction tor in %s', basename)
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms, mask=wmap)

    # if seeing_fwhm is not None:
    #    sex.config['SEEING_FWHM'] = seeing_fwhm * sex.config['PIXEL_SCALE']

    # sex.config['PARAMETERS_LIST'].append('CLASS_STAR')

    # sex.config['CATALOG_NAME'] = 'master-catalogue-i%01d.cat' % step

    LIMIT_AREA = 5000
    idx_small = objects['npix'] < LIMIT_AREA
    objects_small = objects[idx_small]
    NKEEP = 15
    idx_flux = objects_small['flux'].argsort()
    objects_nth = objects_small[idx_flux][-NKEEP:]

    # set of indices of the N first objects
    fluxes = []
    errors = []
    times = []
    airmasses = []

    for idx, frame in enumerate(frames):
        imagename = name_skysub_proc(frame.baselabel, step)

        #sex.config['CATALOG_NAME'] = ('catalogue-%s-i%01d.cat' %
        #                              (frame.baselabel, step))

        # Lauch SExtractor on a FITS file
        # om double image mode
        _logger.info('Runing sextractor in %s', imagename)
        with fits.open(imagename) as hdul:
            header = hdul[0].header
            airmasses.append(header['airmass'])
            times.append(header['tstamp'])
            data_i = hdul[0].data
            data_i = data_i.byteswap().newbyteorder()
            bkg_i = sep.Background(data_i)
            data_sub_i = data_i - bkg_i
        # objects_i = sep.extract(data_sub_i, 1.5, err=bkg_i.globalrms, mask=wmap)
        flux_i, fluxerr_i, flag_i = sep.sum_circle(data_sub_i,
                                                   objects_nth['x'], objects_nth['y'],
                                                   3.0, err=bkg_i.globalrms)


        # Extinction correction
        excor = pow(10, -0.4 * frame.airmass * extinction)

        flux_i = excor * flux_i
        fluxerr_i = excor * fluxerr_i
        fluxes.append(flux_i)
        errors.append(fluxerr_i)

    fluxes_a = numpy.array(fluxes)
    errors_a = numpy.array(errors)
    fluxes_n = fluxes_a / fluxes_a[0]
    errors_a = errors_a / fluxes_a[0]  # sigma
    w = 1.0 / (errors_a) ** 2

    # weighted mean of the flux values
    wdata = numpy.average(fluxes_n, axis=1, weights=w)
    wsigma = 1 / numpy.sqrt(w.sum(axis=1))

    levels = check_photometry_levels
    actions = check_photometry_actions

    x = list(six.moves.range(len(frames)))
    vals, (_, sigma) = check_photometry_categorize(
        x, wdata, levels, tags=actions)
    # n sigma level to plt
    nsig = 3

    if True:
        figure = plt.figure()
        ax = figure.add_subplot(111)

        plot_photometry_check(ax, vals, wsigma, check_photometry_levels, nsig * sigma)
        plt.savefig('figure-relative-flux_i%01d.png' % step)

    for x, _, t in vals:
        try:
            action = _dactions[t]
        except KeyError:
            _logger.warning('Action named %s not recognized, ignoring', t)
            action = default_action
        for p in x:
            action(frames[p])


def check_photometry_categorize(x, y, levels, tags=None):
    '''Put every point in its category.

    levels must be sorted.'''
    x = numpy.asarray(x)
    y = numpy.asarray(y)
    ys = y.copy()
    ys.sort()
    # Mean of the upper half
    m = ys[len(ys) // 2:].mean()
    y /= m
    m = 1.0
    s = ys[len(ys) // 2:].std()
    result = []

    if tags is None:
        tags = list(six.moves.range(len(levels) + 1))

    for l, t in zip(levels, tags):
        indc = y < l
        if indc.any():
            x1 = x[indc]
            y1 = y[indc]
            result.append((x1, y1, t))

            x = x[~indc]
            y = y[~indc]
    else:
        result.append((x, y, tags[-1]))

    return result, (m, s)


def plot_photometry_check(ax, vals, errors, levels, nsigma):
    x = range(len(errors))

    ax.set_title('Relative flux of brightest object')
    for v, c in zip(vals, ['b', 'r', 'g', 'y']):
        ax.scatter(v[0], v[1], c=c)
        w = errors[v[0]]
        ax.errorbar(v[0], v[1], yerr=w, fmt='none', c=c)

    ax.plot([x[0], x[-1]], [1, 1], 'r--')
    ax.plot([x[0], x[-1]], [1 - nsigma, 1 - nsigma], 'b--')
    for f in levels:
        ax.plot([x[0], x[-1]], [f, f], 'g--')

    return ax




def check_position(images_info, sf_data, seeing_fwhm, step=0):
        # FIXME: this method has to be updated

    _logger.info('Checking positions')
    # Check position of bright objects
    weigthmap = 'weights4rms.fits'

    wmap = numpy.zeros_like(sf_data[0])

    # Center of the image
    border = 300
    wmap[border:-border, border:-border] = 1
    fits.writeto(weigthmap, wmap.astype('uint8'), overwrite=True)

    basename = 'result_i%0d.fits' % (step)
    sex = SExtractor()
    sex.config['VERBOSE_TYPE'] = 'QUIET'
    sex.config['PIXEL_SCALE'] = 1
    sex.config['BACK_TYPE'] = 'AUTO'
    if seeing_fwhm is not None and seeing_fwhm > 0:
        sex.config['SEEING_FWHM'] = seeing_fwhm * sex.config['PIXEL_SCALE']
    sex.config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
    sex.config['WEIGHT_IMAGE'] = weigthmap

    sex.config['PARAMETERS_LIST'].append('FLUX_BEST')
    sex.config['PARAMETERS_LIST'].append('FLUXERR_BEST')
    sex.config['PARAMETERS_LIST'].append('FWHM_IMAGE')
    sex.config['PARAMETERS_LIST'].append('CLASS_STAR')

    sex.config['CATALOG_NAME'] = 'master-catalogue-i%01d.cat' % step

    _logger.info('Runing sextractor in %s', basename)
    sex.run('%s,%s' % (basename, basename))

    # Sort catalog by flux
    catalog = sex.catalog()
    catalog = sorted(
        catalog, key=operator.itemgetter('FLUX_BEST'), reverse=True)

    # set of indices of the N first objects
    OBJS_I_KEEP = 10

    # master = [(obj['X_IMAGE'], obj['Y_IMAGE'])
    #    for obj in catalog[:OBJS_I_KEEP]]

    for image in images_info:
        imagename = name_skysub_proc(image.baselabel, step)

        sex.config['CATALOG_NAME'] = ('catalogue-self-%s-i%01d.cat' %
                                      (image.baselabel, step))

        # Lauch SExtractor on a FITS file
        # on double image mode
        _logger.info('Runing sextractor in %s', imagename)
        sex.run(imagename)
        catalog = sex.catalog()

        # data = [(obj['X_IMAGE'], obj['Y_IMAGE']) for obj in catalog]

        # tree = KDTree(data)

        # Search 2 neighbors
        # dists, _ids = tree.query(master, 2, distance_upper_bound=5)

        # for i in dists[:,0]:
        #    print i

        # _logger.info('Mean offset correction for image %s is %f',
        #    imagename, dists[:,0].mean())
        # raw_input('press any key')
