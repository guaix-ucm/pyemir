#
# Copyright 2011-2016 Universidad Complutense de Madrid
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

'''Auxiliary Recipes for EMIR'''

import logging
import math

import numpy
from scipy.ndimage.filters import generic_filter

_logger = logging.getLogger('numina.recipes.emir')


class Trace(object):
    def __init__(self):
        self.tracex = []
        self.tracey = []
        
    def append(self, x, y):
        self.tracex.append(x)
        self.tracey.append(y)
        
    def predict(self, x):
        return self.tracey[-1]

    def reverse(self):
        self.tracex.reverse()
        self.tracey.reverse()


def plot_trace(dd2d, dd1d, deriv, ppmax1, ppmax2, background):
    
#    fig, ax = plt.subplots(2, 2, figsize=(10, 8))
#    ax[0,0].imshow(dd2d)
#    if ppmax1.size > 1:
#        ax[0,0].axhline(ppmax1[0,0], color='green')
#    if ppmax2.size > 1:
#        ax[0,0].axhline(ppmax2[0,0], color='green')
#    ax[0,0].set_axis_off()
#    
#    ax[0,1].plot(dd1d, '*-')
#
#
#
#    ax[1,0].plot(deriv, 'r.-')
#    ax[1,0].axhline(background)
#    ax[1,0].axhline(-background)
#    if ppmax1.size > 1:
#        ax[1,0].plot(ppmax1[:,0], ppmax1[:,1], 'og')
#    if ppmax2.size > 1:
#        ax[1,0].plot(ppmax2[:,0], -ppmax2[:,1], 'ob');
#    plt.show()
#        #plt.close(fig)
    pass


WW = numpy.array([[0.5, -1.0, 0.5], [-0.5, 0.0, 0.5], [0.0, 1.0, 0.0]])


def peak_detection_basic(y, x=None, background=0.0):
    '''Basic peak detection.'''

    peaks = []
    
    for idx, value in enumerate(y[1:-1], start=1):
        if (value > y[idx-1]) and (value > y[idx+1]) and (value > background):
            peaks.append((idx, value))
    
    
    return numpy.array(peaks)


def interp_3(dd):
    
    coeff= numpy.dot(WW, dd)
    
    return coeff


def wc_to_pix(x):
    return int(math.floor(x + 0.5))


def img_norm(img):
    '''Normalize image between +1 and -1.'''
    mn = img.min()
    mx = img.max()
    if mn == mx:
        return numpy.zeros_like(img)

    res = 2 * (img - mn) / (mx - mn) - 1
    return res


def _internal_trace(img, trace1, trace2, col, step, hs, ws, tol=2, direction=1, doplot=False):

    tolcounter = tol
    axis_size = img.shape[1]
    
    # first derivatice cuadratic fit
    _w1 = [-2.0, -1.0, 0.0, 1.0, 2.0]
    _w2 = [-0.1071, -0.0714, -0.0357, 0.0, 0.0357, 0.0714, 0.1071]
    
    while (col+ step > hs) and (col + step +hs) < axis_size:
        #print 'we are in column', col
        col += direction * step
        #print 'we go to column', col
        prediction1 = trace1.predict(col)
        prediction2 = trace2.predict(col)
        #print 'we predict the borders will be in coordinates', prediction1, prediction2
        pred_pix1 = wc_to_pix(prediction1)
        pred_pix2 = wc_to_pix(prediction2)
    
        reg1 = pred_pix1 - ws
        reg2 = pred_pix2 + ws
        #print 'region is', reg1, reg2, col-hs,col+hs

        dd2d = img[reg1:reg2+1, col-hs:col+hs +1]
    
        # Colapse
        dd1d = dd2d.mean(axis=1)
        # Derivative
        
        def filter_fun(a):
            return numpy.sum(a * _w2)
        
        out = generic_filter(dd1d, filter_fun, size=len(_w2))
        # Quantiles to estimate background noise
        q75, _q50, q25 = numpy.percentile(out, [75, 50, 25])
        iqr = q75 - q25
        background = 2.5 * iqr

        # peaks, positive (first border)
        ppmax1 = peak_detection_basic(out, background=background)
        # peaks, negative (second border)
        ppmax2 = peak_detection_basic(-out, background=background)
    

        # find nearest peak to prediction
    
        # check the peak is not further than npixels'

        # Do consistency check in the peak position (left positive, right negative, etc)
    
        if doplot:
            plot_trace(dd2d, dd1d, out, ppmax1, ppmax2, background)
        
        if ppmax1.size < 1 or ppmax2.size < 1:
            if tolcounter > 1:
                #print 'border missing, try again'
                tolcounter -= 1
                continue
            else:
                #print "border missing, no more tries"
                return trace1, trace2
    
        tolcounter = tol

        p1 = ppmax1[0,0]
        p2 = ppmax2[0,0]
    
        # fit the peak to a parabole
        # fine tuning
        A, B, C = interp_3(out[p1-1:p1+2])
        e1, _e2 = -B / (2*A), C - B*B / (4*A)
        trace1.append(col, reg1 + p1 + e1)

        A, B, C = interp_3(-out[p2-1:p2+2])
        e1, _e2 = -B / (2*A), C - B*B / (4*A)    
        trace2.append(col, reg1 + p2 + e1)
        
    return trace1, trace2


def _internal_trace_c(img, trace1, trace2, col, step, hs, ws, tol, doplot=False):

    _internal_trace(img, trace1, trace2, col, step, hs, ws, tol, direction=-1, doplot=doplot)

    trace1.reverse()
    trace2.reverse()
    
    _internal_trace(img, trace1, trace2, col, step, hs, ws, tol, direction=+1, doplot=doplot)
    
    return trace1, trace2


def ex_region(img, x, y1, y2, step, hs, ws, tol=2, doplot=False):
    
    trace1 = Trace()
    trace2 = Trace()
    # Initial point
    trace1.append(x, y1)
    trace2.append(x, y2)
    
    _internal_trace_c(img, trace1, trace2, x, step, hs, ws, tol, doplot=doplot)
    
    pfit1 = numpy.polyfit(trace1.tracex, trace1.tracey, deg=5)
    pfit2 = numpy.polyfit(trace2.tracex, trace2.tracey, deg=5)
    
    return (trace1.tracex[0],
            trace1.tracex[-1],
            wc_to_pix(numpy.min(trace1.tracey)),
            wc_to_pix(numpy.max(trace2.tracey)),
            pfit1,
            pfit2
           )


def convert_to_(x, y, ax, ay):
    col = wc_to_pix(x - 1)
    y1 = y - 1 - 0.5 * ay
    y2 = y - 1 + 0.5 * ay
    return col, y1, y2
