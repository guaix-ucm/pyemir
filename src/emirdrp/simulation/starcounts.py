#
# Copyright 2008-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import math

import importlib_resources

import numpy
import scipy.interpolate as sil


class Counts:
    INTEGRAL_COUNTS = (0,)
    DIFFERENTIAL_COUNTS = 1


class PhotometricFilter:
    FILTER_NONE = 0
    FILTER_B = 1
    FILTER_V = 2
    FILTER_Z = 3
    FILTER_J = 4
    FILTER_H = 5
    FILTER_K = 6


class BaseModel:
    def __init__(self, name):
        self.name = name


class RBModel(BaseModel):
    """Star counts model from Ratnatunga & Bahcall.

    Number of stars per square arc minute
    """

    def __init__(self):
        super().__init__("Ratnatunga & Bahcall")

        _counts = numpy.array(
            [
                [12, 14, 16, 18, 20, 22, 24, 26, 28],
                [0.0, 1.3e-2, 1.8e-2, 4.8e-2, 1.1e-1, 9.0e-2, 2.7e-2, 3.4e-3, 1.6e-4],
                [0.0, 8.2e-3, 2.6e-2, 2.6e-2, 4.5e-2, 1.2e-1, 1.1e-1, 3.9e-2, 5.6e-3],
                [0.0, 2.0e-3, 1.0e-2, 4.4e-2, 1.1e-1, 2.2e-1, 4.6e-1, 8.1e-1, 1.1e-0],
            ]
        )
        # Order 3 spline
        self._spl = [sil.splrep(_counts[0], _counts[i]) for i in (1, 2, 3)]
        self._min = _counts[0, 0]

        self._null = [0.0, 0.0, 0.0]
        self._colorsVK = [4.7, 2.5, 1.4]

        self.filters = [PhotometricFilter.FILTER_K]

    def differential_counts(self, mag, filter=PhotometricFilter.FILTER_K):
        if filter == PhotometricFilter.FILTER_V:
            return self._differential_counts_color(mag, self._null)
        elif filter == PhotometricFilter.FILTER_K:
            return self._differential_counts_color(mag, self._colorsVK)

    def integral_counts(self, mag, filter=PhotometricFilter.FILTER_K):
        if filter == PhotometricFilter.FILTER_K:
            return self._integral_counts_color(mag, self._colorsVK)
        return 0.0

    def _integral_counts_color(self, mag, colors):
        def subint(i):
            return max(
                sil.splint(self._min + colors[i], mag + colors[i], self._spl[i]), 0.0
            )

        return sum(subint(i) for i in (0, 1, 2))

    def _differential_counts_color(self, mag, colors):
        def subspl(i):
            return max(sil.splev(mag, self._spl[i]), 0.0)

        return sum(subspl(i) for i in (0, 1, 2))


class SpagnaModel(BaseModel):
    """Star counts model from Spagna 1999.

    Number of stars per square arc minute
    """

    def __init__(self):
        super().__init__("Spagna 1999")

        res = importlib_resources.files("emirdrp.simulation").joinpath("spagna-J.dat")
        with open(res) as fd:
            _J_counts_data = numpy.loadtxt(fd)

            # Data in file is for square degree
            _J_counts_data[:, 1:3] /= 3600.0
            self._spl_J_1 = sil.splrep(_J_counts_data[:, 0], _J_counts_data[:, 1])
            self._spl_J_2 = sil.splrep(_J_counts_data[:, 0], _J_counts_data[:, 2])

        res = importlib_resources.files("emirdrp.simulation").joinpath("spagna-K.dat")
        with open(res) as fd:
            _K_counts_data = numpy.loadtxt(fd)
            # Data in file is for square degree
            _K_counts_data[:, 1:3] /= 3600.0
            self._spl_K_1 = sil.splrep(_K_counts_data[:, 0], _K_counts_data[:, 1])
            self._spl_K_2 = sil.splrep(_K_counts_data[:, 0], _K_counts_data[:, 2])

        self.filters = [PhotometricFilter.FILTER_K, PhotometricFilter.FILTER_J]

    def integral_counts(self, mag, filter=PhotometricFilter.FILTER_K):
        if filter == PhotometricFilter.FILTER_J:
            return sil.splev(mag, self._spl_J_2)
        elif filter == PhotometricFilter.FILTER_K:
            return sil.splev(mag, self._spl_K_2)
        else:
            return 0.0

    def differential_counts(self, mag, filter=PhotometricFilter.FILTER_K):
        if filter == PhotometricFilter.FILTER_J:
            return sil.splev(mag, self._spl_J_1) / 3600.0
        elif filter == PhotometricFilter.FILTER_K:
            return sil.splev(mag, self._spl_K_1) / 3600.0
        else:
            return 0.0


class BSModel(BaseModel):
    """Star counts model from Bahcall & Soneira.

    Number of stars per square arc minute
    """

    def __init__(self):
        super().__init__("Bahcall & Soneira")
        self.params = numpy.array(
            [
                [200, -0.2, 0.01, 2, 15, 400, -0.26, 0.065, 1.5, 17.5],
                [925, -0.132, 0.035, 3, 15.75, 1050, -0.18, 0.087, 2.5, 17.5],
                [235, -0.227, 0.0, 1.5, 17, 370, -0.175, 0.06, 2.0, 18],
                [950, -0.124, 0.027, 3.1, 16.60, 910, -0.167, 0.083, 2.5, 18],
            ]
        )

    def integral_counts(self, mag, filter=PhotometricFilter.FILTER_V):
        if filter == PhotometricFilter.FILTER_V:
            return self.dFunction(mag, 0.0, 0.5 * math.pi, *self.params[1])
        elif filter == PhotometricFilter.FILTER_B:
            return self.dFunction(mag, 0.0, 0.5 * math.pi, *self.params[4])
        else:
            return 0.0

    def differential_counts(self, mag, filter=PhotometricFilter.FILTER_V):
        if filter == PhotometricFilter.FILTER_V:
            return self.dFunction(mag, 0.0, 0.5 * math.pi, *self.params[0])
        elif filter == PhotometricFilter.FILTER_B:
            return self.dFunction(mag, 0.0, 0.5 * math.pi, *self.params[3])
        else:
            return 0.0

    def mu(self, magnitude):
        """Compute parameter \\mu of Appendix B."""
        if magnitude <= 12.0:
            return 0.03
        if magnitude <= 20:
            return 0.0075 * (magnitude - 12.0) + 0.03
        # if magnitude > 20:
        return 0.09

    def gamma(self, magnitude):
        if magnitude <= 12.0:
            return 0.36
        if magnitude <= 20:
            return 0.04 * (12.0 - magnitude) + 0.36
        # if magnitude > 20:
        return 0.04

    def sigma(self, galacticLongitude, galacticLatitude):
        return 1.45 - 0.20 * math.cos(galacticLongitude) * math.cos(galacticLatitude)

    def dFunction(
        self,
        mag,
        glongitude,
        glatitude,
        C1,
        alpha,
        beta,
        delta,
        mStar,
        C2,
        kappa,
        eta,
        lambdax,
        mDagger,
    ):
        firstTerm = (
            C1
            * pow(10.0, beta * (mag - mStar))
            / pow(1 + pow(10.0, alpha * (mag - mStar)), delta)
            / pow(
                math.sin(glatitude)
                * (1 - self.mu(mag) / math.tan(glatitude) * math.cos(glongitude)),
                3.0 - 5 * self.gamma(mag),
            )
        )
        secondTerm = (
            C2
            * pow(10.0, eta * (mag - mDagger))
            / pow(1 + pow(10.0, kappa * (mag - mDagger)), lambdax)
            / pow(
                1 - math.cos(glatitude) * math.cos(glongitude),
                self.sigma(glongitude, glatitude),
            )
        )
        return (firstTerm + secondTerm) / 3600.0


if __name__ == "__main__":
    import numpy as np
    import numpy.random
    from math import sqrt, log

    rbmodel = RBModel()
    # print rbmodel.integral_counts(18)
    # print rbmodel.differential_counts(18)
    sgmodel = SpagnaModel()
    # print sgmodel.integral_counts(18)
    # print sgmodel.differential_counts(18)
    bsmodel = BSModel()
    # print sgmodel.integral_counts(18)
    # print bsmodel.differential_counts(18)

    class GeneralRandom:
        """Recipe from http://code.activestate.com/recipes/576556/."""

        def __init__(self, x, p, Nrl=1000):
            self.x = x
            self.pdf = p / p.sum()
            self.cdf = self.pdf.cumsum()
            self.inversecdfbins = Nrl
            self.Nrl = Nrl
            y = np.arange(Nrl) / float(Nrl)
            # _delta = 1.0 / Nrl
            self.inversecdf = np.zeros(Nrl)
            self.inversecdf[0] = self.x[0]
            cdf_idx = 0
            for n in range(1, self.inversecdfbins):
                while self.cdf[cdf_idx] < y[n] and cdf_idx < Nrl:
                    cdf_idx += 1
                # Seems a linear interpolation
                self.inversecdf[n] = self.x[cdf_idx - 1] + (
                    self.x[cdf_idx] - self.x[cdf_idx - 1]
                ) * (y[n] - self.cdf[cdf_idx - 1]) / (
                    self.cdf[cdf_idx] - self.cdf[cdf_idx - 1]
                )
                if cdf_idx >= Nrl:
                    break
            self.delta_inversecdf = np.concatenate((np.diff(self.inversecdf), [0]))

        def random(self, N=1000):
            idx_f = np.random.uniform(size=N, high=self.Nrl - 1)
            idx = np.array(idx_f, "i")
            y = self.inversecdf[idx] + (idx_f - idx) * self.delta_inversecdf[idx]
            return y

    scmodel = sgmodel
    magmin = 12.0
    magmax = 25.0

    plate_scale = 0.2
    detector_shape = (2048, 2048)
    pixel_area = detector_shape
    detector_area = (plate_scale**2 * detector_shape[0] * detector_shape[1]) / 3600.0

    nstars = int(
        round(
            detector_area
            * (scmodel.integral_counts(magmax) - scmodel.integral_counts(magmin))
        )
    )
    print("nstars", nstars)

    seeing = 1.0

    scale = 2 * sqrt(2 * log(2))

    delta = 0.1
    steps = int(round((magmax - magmin) / delta))
    x = magmin + delta * np.arange(steps)
    p = np.array([scmodel.differential_counts(i) for i in x])

    g = GeneralRandom(x, p)

    np.random.seed(100)
    y = np.random.uniform(high=pixel_area[0], size=nstars)
    x = np.random.uniform(high=pixel_area[1], size=nstars)
    sigmas = [seeing / scale / plate_scale] * nstars
    ints = [100] * nstars
    mag = g.random(N=nstars)
