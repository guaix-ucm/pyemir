# Version 29 May 2015
#------------------------------------------------------------------------------

from __future__ import division
from __future__ import print_function

import sys

import six
from six.moves import input
import numpy as np
from numpy.polynomial import polynomial
import matplotlib.pyplot as plt

from .zscale import zscale

#------------------------------------------------------------------------------

class Slitlet(object):
    """Slitlet definition."""

    def __init__(self,bb_nc1,bb_nc2,bb_ns1,bb_ns2):
        """Initialize slitlet with bounding box.

        Note that bb_nc1, bb_nc2, bb_ns1 and bb_ns2 are in pixel units.
        """
        # protections
        if type(bb_nc1) not in [np.int, np.int64]:
            raise ValueError('bb_nc1='+str(bb_nc1)+' must be integer')
        if type(bb_nc2) not in [np.int, np.int64]:
            raise ValueError('bb_nc2='+str(bb_nc2)+' must be integer')
        if type(bb_ns1) not in [np.int, np.int64]:
            raise ValueError('bb_ns1='+str(bb_ns1)+' must be integer')
        if type(bb_ns2) not in [np.int, np.int64]:
            raise ValueError('bb_ns2='+str(bb_ns2)+' must be integer')
        # define enclosing bounding box (in pixel units)
        self.nc_lower_bb_pix = bb_nc1
        self.nc_upper_bb_pix = bb_nc2
        self.ns_lower_bb_pix = bb_ns1
        self.ns_upper_bb_pix = bb_ns2
        # define polynomials for boundaries (in pixel units)
        self.set_nc_coeff_lower_boundary_pix(self.nc_lower_bb_pix)
        self.set_nc_coeff_upper_boundary_pix(self.nc_upper_bb_pix)
        self.set_ns_coeff_lower_boundary_pix(self.ns_lower_bb_pix)
        self.set_ns_coeff_upper_boundary_pix(self.ns_upper_bb_pix)

    def set_nc_coeff_lower_boundary_pix(self, coeffs):
        self.nc_coeff_lower_boundary_pix = np.array(coeffs)
        self.nc_funct_lower_boundary_pix = \
          polynomial.Polynomial(self.nc_coeff_lower_boundary_pix)

    def set_nc_coeff_upper_boundary_pix(self, coeffs):
        self.nc_coeff_upper_boundary_pix = np.array(coeffs)
        self.nc_funct_upper_boundary_pix = \
          polynomial.Polynomial(self.nc_coeff_upper_boundary_pix)

    def set_ns_coeff_lower_boundary_pix(self, coeffs):
        self.ns_coeff_lower_boundary_pix = np.array(coeffs)
        self.ns_funct_lower_boundary_pix = \
          polynomial.Polynomial(self.ns_coeff_lower_boundary_pix)

    def set_ns_coeff_upper_boundary_pix(self, coeffs):
        self.ns_coeff_upper_boundary_pix = np.array(coeffs)
        self.ns_funct_upper_boundary_pix = \
          polynomial.Polynomial(self.ns_coeff_upper_boundary_pix)

    def extract_midsp(self, image2d, nwidth, method='sum', LDEBUG=False):
        """Extract a spectrum from the central track within a fixed window.

        The central track is computed as the mean location between the lower
        and the upper boundaries. The function extracts the central spectrum by
        coadding all the flux in a vertical window of 'nwidth' scans, using
        fractions of pixel at the borders when necessary.

        Parameters
        ----------
        image2d : 2d numpy array (float)
            Input image.
        nwidth : int
            Number of scans to be coadded.
        method : string
            Method to be used to compute spectrum:
            - sum: sum all the flux
            - mean: sum all the flux and divide by the total number of coadded
              scans
        LDEBUG : bool
            If True the function plots and prints out additional information.

        Returns
        -------
        xsp_float : 1d numpy array (float)
            X-coordinate (channel number) of the pixels of the extracted
            spectrum.
        fluxsp : 1d numpy array (float)
            Flux of the extracted spectrum.
        """

        # protections
        if method not in ['sum', 'mean']:
            raise ValueError('Invalid method='+method)
        if type(nwidth) not in [np.int, np.int64]:
            raise ValueError('nwidth='+str(nwidth)+' must be integer')
        if nwidth < 1:
            raise ValueError('nwidth must be >= 1')

        # bounding box
        nc1 = self.nc_lower_bb_pix # int
        nc2 = self.nc_upper_bb_pix # int
        ns1 = self.ns_lower_bb_pix # int
        ns2 = self.ns_upper_bb_pix # int

        # array with channel values
        xsp = np.array(range(nc1,nc2+1)) # int

        # evaluate lower and upper boundaries in each channel
        ybound1 = self.ns_funct_lower_boundary_pix(xsp) # float
        ybound2 = self.ns_funct_upper_boundary_pix(xsp) # float

        # compute central track
        yspcen = (ybound1+ybound2)/2 # float

        # compute upper and lower limits around the central track
        yspmin = yspcen-float(nwidth)/2 # float
        yspmax = yspcen+float(nwidth)/2 # float

        # clip previous limits if outside bounding box 
        # (in-place clip using the parameter out=...)
        np.clip(yspmin, float(ns1)-0.5, float(ns2)+0.5, out=yspmin)
        np.clip(yspmax, float(ns1)-0.5, float(ns2)+0.5, out=yspmax)

        # coadd all pixels between previous limits, using fractions of pixels
        # when necessary
        nchaneff = nc2-nc1+1
        fluxsp = np.zeros(nchaneff, dtype=np.float)
        for j in range(nchaneff):
            fluxsp[j] = 0.0
            jdum = xsp[j]-1  # array coordinate (channel_number - 1), int
            fdum1 = yspmin[j]-1 # array coordinate (scan_number - 1), float
            fdum2 = yspmax[j]-1 # array coordinate (scan_number - 1), float
            if fdum2 <= fdum1:
                pass # nothing to sum
            else:
                idum1 = int(fdum1+0.5) # array coordinate, int
                idum2 = int(fdum2+0.5) # array coordinate, int
                if idum2 == idum1: # sum fraction of the same scan
                    fluxsp[j] += image2d[idum1,jdum]*(fdum2-fdum1)
                else:
                    fraction1 = float(idum1)-fdum1+0.5 # fraction first scan
                    fluxsp[j] += image2d[idum1,jdum]*fraction1
                    fraction2 = fdum2-float(idum2)+0.5 # fraction last scan
                    fluxsp[j] += image2d[idum2,jdum]*fraction2
                    if idum2 > idum1+1: # full scans in between
                        fluxsp[j] += np.sum(image2d[idum1+1:idum2,jdum],0)

        # if requested, obtain average per scan
        if method == 'mean':
            deltaysp = yspmax-yspmin
            fluxsp /= deltaysp

        if LDEBUG:
            # plot image with bounding box, boundaries, etc.
            naxis2, naxis1 = image2d.shape
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.autoscale(False)
            ax.set_xlim([0.5,float(naxis1)+0.5])
            ax.set_ylim([0.5,float(naxis2)+0.5])
            # plot image
            z1, z2 = zscale(image2d)
            im = plt.imshow(image2d, cmap='hot', aspect='auto',
                            vmin = z1, vmax = z2,
                            interpolation = 'nearest', origin = 'low',
                            extent = [0.5,float(naxis1)+0.5,
                                      0.5,float(naxis2)+0.5])
            # plot bounding box
            xplot = [nc1,nc1,nc2,nc2,nc1]
            yplot = [ns1,ns2,ns2,ns1,ns1]
            ax.plot(xplot, yplot, 'y-')
            # plot boundaries
            ax.plot(xsp,ybound1,'m-')
            ax.plot(xsp,ybound2,'m-')
            ax.plot(xsp,yspcen,'g-')
            ax.plot(xsp,yspmin,'c-')
            ax.plot(xsp,yspmax,'c-')
            # show plot
            plt.show(block=False)
            input('press <RETURN> to continue...')

            # plot extracted spectrum
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlim([1,naxis1])
            ax.plot(xsp,fluxsp,'k-')
            plt.show(block=False)
            input('press <RETURN> to continue...')

        # return result
        xsp_float = xsp.astype(np.float)
        return xsp_float, fluxsp

