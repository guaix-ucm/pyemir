#
# Copyright 2011-2012 Universidad Complutense de Madrid
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

'''

Routines shared by image mode recipes

'''

import logging
import operator

import numpy
import pyfits
from numina.util.sextractor import SExtractor

from .naming import name_skysub_proc

_logger = logging.getLogger('numina.recipes.emir')

# Actions to carry over images when checking the flux
# of the objects in different images
def warn_action(img):
    _logger.warn('Image %s has low flux in objects', img.baselabel)
    img.valid_science = True

def reject_action(img):
    img.valid_science = False
    _logger.info('Image %s rejected, has low flux in objects', img.baselabel)            
    pass

def default_action(img):
    _logger.info('Image %s accepted, has correct flux in objects', img.baselabel)      
    img.valid_science = True

# Actions
_dactions = {'warn': warn_action, 'reject': reject_action, 'default': default_action}

def check_photometry(frames, sf_data, seeing_fwhm, step=0, border=300, extinction=0.0,
                     check_photometry_levels=[0.5, 0.8], check_photometry_actions=['warn', 'warn', 'default'], figure=None):
    # Check photometry of few objects
    weigthmap = 'weights4rms.fits'
    
    wmap = numpy.zeros_like(sf_data[0])
    
    # Center of the image
    wmap[border:-border, border:-border] = 1                    
    pyfits.writeto(weigthmap, wmap, clobber=True)
    
    basename = 'result_i%0d.fits' % (step)
    sex = SExtractor()
    sex.config['VERBOSE_TYPE'] = 'QUIET'
    sex.config['PIXEL_SCALE'] = 1
    sex.config['BACK_TYPE'] = 'AUTO' 
    if seeing_fwhm is not None:
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
    catalog = sorted(catalog, key=operator.itemgetter('FLUX_BEST'), reverse=True)
    
    # set of indices of the N first objects
    OBJS_I_KEEP = 3
    indices = set(obj['NUMBER'] for obj in catalog[:OBJS_I_KEEP])
    
    base = numpy.empty((len(frames), OBJS_I_KEEP))
    error = numpy.empty((len(frames), OBJS_I_KEEP))
    
    for idx, frame in enumerate(frames):
        imagename = name_skysub_proc(frame.baselabel, step)

        sex.config['CATALOG_NAME'] = 'catalogue-%s-i%01d.cat' % (frame.baselabel, step)

        # Lauch SExtractor on a FITS file
        # om double image mode
        _logger.info('Runing sextractor in %s', imagename)
        sex.run('%s,%s' % (basename, imagename))
        catalog = sex.catalog()
        
        # Extinction correction
        excor = pow(10, -0.4 * frame.airmass * extinction)
        base[idx] = [obj['FLUX_BEST'] / excor
                                 for obj in catalog if obj['NUMBER'] in indices]
        error[idx] = [obj['FLUXERR_BEST'] / excor
                                 for obj in catalog if obj['NUMBER'] in indices]
    
    data = base / base[0]
    err = error / base[0] # sigma
    w = 1 / err / err
    # weighted mean of the flux values
    wdata = numpy.average(data, axis=1, weights=w)
    wsigma = 1 / numpy.sqrt(w.sum(axis=1))
    
    levels = check_photometry_levels
    actions = check_photometry_actions
    
    x = range(len(frames))
    vals, (_, sigma) = check_photometry_categorize(x, wdata, levels, tags=actions)
    # n sigma level to plt
    n = 3
    
    if figure is not None:
        check_photometry_plot(figure, vals, wsigma, levels, n * sigma)
    
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
    m = ys[len(ys) / 2:].mean()
    y /= m
    m = 1.0
    s = ys[len(ys) / 2:].std()
    result = []

    if tags is None:
        tags = range(len(levels) + 1)

    for l, t in zip(levels, tags):
        indc = y < l
        if indc.any():
            x1 = x[indc]
            y1 = y[indc]
            result.append((x1, y1, t))

            x = x[indc == False]
            y = y[indc == False]
    else:
        result.append((x, y, tags[-1]))

    return result, (m, s)
           
def check_photometry_plot(figure, vals, errors, levels, nsigma, step=0):
    x = range(len(errors))
    figure.clf()
    ax = figure.add_subplot(111)
    ax.set_title('Relative flux of brightest object')
    for v,c in zip(vals, ['b', 'r', 'g', 'y']):
        ax.scatter(v[0], v[1], c=c)
        w = errors[v[0]]
        ax.errorbar(v[0], v[1], yerr=w, fmt=None, c=c)
        

    ax.plot([x[0], x[-1]], [1, 1], 'r--')
    ax.plot([x[0], x[-1]], [1 - nsigma, 1 - nsigma], 'b--')
    for f in levels:
        ax.plot([x[0], x[-1]], [f, f], 'g--')
        
    figure.canvas.draw()
    figure.savefig('figure-relative-flux_i%01d.png' % step)

def check_position(images_info, sf_data, seeing_fwhm, step=0):
        # FIXME: this method has to be updated

        _logger.info('Checking positions')
        # Check position of bright objects
        weigthmap = 'weights4rms.fits'
        
        wmap = numpy.zeros_like(sf_data[0])
        
        # Center of the image
        border = 300
        wmap[border:-border, border:-border] = 1                    
        pyfits.writeto(weigthmap, wmap, clobber=True)
        
        basename = 'result_i%0d.fits' % (step)
        sex = SExtractor()
        sex.config['VERBOSE_TYPE'] = 'QUIET'
        sex.config['PIXEL_SCALE'] = 1
        sex.config['BACK_TYPE'] = 'AUTO'
        if  seeing_fwhm is not None and seeing_fwhm > 0:
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
        catalog = sorted(catalog, key=operator.itemgetter('FLUX_BEST'), reverse=True)
        
        # set of indices of the N first objects
        OBJS_I_KEEP = 10
        
        # master = [(obj['X_IMAGE'], obj['Y_IMAGE']) for obj in catalog[:OBJS_I_KEEP]]
        
        for image in images_info:
            imagename = name_skysub_proc(image.baselabel, step)

            sex.config['CATALOG_NAME'] = 'catalogue-self-%s-i%01d.cat' % (image.baselabel, step)

            # Lauch SExtractor on a FITS file
            # on double image mode
            _logger.info('Runing sextractor in %s', imagename)
            sex.run(imagename)
            catalog = sex.catalog()
            
            
            # data = [(obj['X_IMAGE'], obj['Y_IMAGE']) for obj in catalog]
            
            #tree = KDTree(data)
            
            # Search 2 neighbors
            #dists, _ids = tree.query(master, 2, distance_upper_bound=5)
            
            #for i in dists[:,0]:
            #    print i
            
            
            #_logger.info('Mean offset correction for image %s is %f', imagename, dists[:,0].mean())
            #raw_input('press any key')
