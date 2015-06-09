# Version 29 May 2015
#------------------------------------------------------------------------------

from __future__ import division
from __future__ import print_function

import six
from six.moves import input

import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------------

def findPeaks_spectrum(sx, nwinwidth, data_threshold = 0, 
                       LDEBUG = False, LPLOT = False):
    """Find peaks in the 1d-numpy array sx.

    Note that nwinwidth must be an odd number. The function imposes that the
    fluxes in the nwinwidth/2 points to the left (and right) of the peak
    decrease monotonously as one moves away from the peak.

    Parameters
    ----------
    sx : 1d numpy array (float)
        Input 1D spectrum.
    nwinwidth : int
        Width of the window where the peak must be found. This number must be
        odd.
    data_threshold : float
        Minimum signal in the peak (optional).
    LDEBUG : bool
        If True the function prints out additional information.
    LPLOT : bool
        If True the function plots the spectrum and the peaks.

    Returns
    -------
    ipeaks : 1d numpy array (int)
        Indices of the input array sx in which the peaks have been found. Note
        that these numbers are not channels, but array indices starting from
        zero.
        
    """

    sx_shape = sx.shape
    nmed = nwinwidth//2

    if LDEBUG:
        print('findPeaks_spectrum> sx shape......:',sx_shape)
        print('findPeaks_spectrum> nwinwidth.....:',nwinwidth)
        print('findPeaks_spectrum> nmed..........:',nmed)
        print('findPeaks_spectrum> data_threshold:',data_threshold)
        print('findPeaks_spectrum> the first and last', nmed, \
              'pixels will be ignored')

    ipeaks = []

    if sx_shape[0] < nwinwidth:
        raise ValueError('findPeaks_spectrum> ERROR: invalid nwinwidth')
        return np.array(ipeaks)

    i = nmed
    while i < sx_shape[0]-nmed:
        if sx[i] > data_threshold:
            lpeakgood = True

            j = 0
            loop = True
            while loop:
                if sx[i-nmed+j] > sx[i-nmed+j+1]: lpeakgood = False
                j += 1
                loop = (j < nmed) and lpeakgood

            if lpeakgood:
                j = nmed+1
                loop = True
                while loop:
                    if sx[i-nmed+j-1] < sx[i-nmed+j]: lpeakgood = False
                    j += 1
                    loop = (j < nwinwidth) and lpeakgood
                
            if lpeakgood: 
                ipeaks.append(i)
                i += nwinwidth-1
            else:
                i += 1
        else:
            i += 1

    if LDEBUG:
        npeaks = len(ipeaks)
        print('findPeaks_spectrum> number of peaks found:',npeaks)
        for val in zip(range(npeaks), ipeaks):
            print(val)

    if LPLOT:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(list(six.moves.range(sx_shape[0])), sx, 'k-')
        ax.set_xlabel('array index')
        ax.set_ylabel('spectrum intensity')
        ax.set_title('function findPeaks_spectrum')
        plt.show(block=False)
        input('press <RETURN> to continue...')

    return np.array(ipeaks, dtype=np.int)

#------------------------------------------------------------------------------

def refinePeaks_spectrum(sx, ipeaks, nwinwidth, method=2, LDEBUG=False):
    """Refine the peak location previously found by findPeaks_spectrum()

    Parameters
    ----------
    sx : 1d numpy array, float
        Input 1D spectrum.
    ipeaks : 1d numpy array (int)
        Indices of the input array sx in which the peaks were initially found.
    nwinwidth : int
        Width of the window where the peak must be found.
    method : int
        Indicates which function will be employed to refine the peak location.
        method  = 1 -> fit to 2nd order polynomial
        method != 1 -> fit to gaussian
    LDEBUG : bool
        If True the function plots and prints out additional information.

    Returns
    -------
    xfpeaks : 1d numpy array (float)
        X-coordinates in which the refined peaks have been found.
 
    """
    nmed = nwinwidth//2

    xfpeaks = np.zeros(len(ipeaks))

    for iline in range(len(ipeaks)):
        jmax = ipeaks[iline]
        sx_peak_flux = sx[jmax]
        x_fit = np.zeros(0)
        y_fit = np.zeros(0)
        for j in range(jmax-nmed,jmax+nmed+1):
            x_fit = np.concatenate((x_fit, np.array([float(j-jmax)])))
            y_fit = np.concatenate((y_fit, np.array([sx[j]/sx_peak_flux])))

        if method == 1:
            poly = np.polyfit(x_fit, y_fit, 2)
            refined_peak = -poly[1]/(2.0*poly[0])+jmax
        else:
            poly = np.polyfit(x_fit, np.log(y_fit), 2)
            A = np.exp(poly[2]-poly[1]*poly[1]/(4*poly[0]))
            x0 = -poly[1]/(2*poly[0])
            sigma = np.sqrt(-1/(2*poly[0]))
            refined_peak = x0+jmax

        xfpeaks[iline] = refined_peak

        if LDEBUG:
            plt.figure()
            xmin = x_fit.min()-1
            xmax = x_fit.max()+1
            ymin = 0
            ymax = y_fit.max()*1.10
            plt.axis([xmin,xmax,ymin,ymax])
            plt.xlabel('channel (around initial integer peak)')
            plt.ylabel('Normalized no. of counts')
            plt.title('Fit to line at channel '+str(jmax))
            plt.plot(x_fit,y_fit,"bo")
            n_plot = 1000
            x_plot = np.zeros(n_plot)
            y_plot = np.zeros(n_plot)
            for j in range(n_plot):
                x_plot[j] = -nmed+float(j)/float(n_plot-1)*2*nmed
                if method == 1:
                    y_plot[j] = poly[2]+ \
                                poly[1]*x_plot[j]+ \
                                poly[0]*x_plot[j]*x_plot[j]
                else:
                    y_plot[j] = A*np.exp(-(x_plot[j]-x0)**2/(2*sigma**2))
            plt.plot(x_plot,y_plot,color="red")
            plt.show(block=False)
            print('refined_peak:',refined_peak)
            answer = input('Press <CR> to continue...')
            plt.close()

    return xfpeaks

