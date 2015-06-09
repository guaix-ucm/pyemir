# Version 29 May 2015
#------------------------------------------------------------------------------

from __future__ import division
from __future__ import print_function

import numpy as np

#------------------------------------------------------------------------------

def sigmaG(x, LDEBUG = False):
    """Compute a robust estimator of the standard deviation

    See Eq. 3.36 (page 84) in Statistics, Data Mining, and Machine
    in Astronomy, by Ivezic, Connolly, VanderPlas & Gray

    Parameters
    ----------
    x : 1d numpy array, float
        Array of input values which standard deviation is requested.
    LDEBUG : bool
        If True prints values.

    Returns
    -------
    sigmag : float
        Robust estimator of the standar deviation
    """

    q25, q75 = np.percentile(x,[25.0, 75.0])
    sigmag = 0.7413*(q75-q25)

    if LDEBUG:
        print('>>> Robust standard deviation:',sigmag)

    return sigmag


#------------------------------------------------------------------------------

def statsummary(x, LDEBUG = False):
    """Compute basic statistical parameters

    Parameters
    ----------
    x : 1d numpy array, float
        Input array with values which statistical properties are requested.
    LDEBUG : bool
        If True prints values.

    Returns
    -------
    result : tuple, floats
        Minimum, percentile 25, percentile 50, mean, percentile 75, maximum,
        standard deviation, and robust standard deviation.

    """

    result = np.min(x),                 \
             np.percentile(x,25),       \
             np.percentile(x,50),       \
             np.mean(x),                \
             np.percentile(x,75),       \
             np.max(x),                 \
             np.std(x),                 \
             sigmaG(x),                 \
             np.percentile(x,15.86553), \
             np.percentile(x,84.13447)

    if LDEBUG:
        print('>>> Minimum..................:',result[0])
        print('>>> 1st Quartile.............:',result[1])
        print('>>> Median...................:',result[2])
        print('>>> Mean.....................:',result[3])
        print('>>> 3rd Quartile.............:',result[4])
        print('>>> Maximum..................:',result[5])
        print('>>> Standard deviation.......:',result[6])
        print('>>> Robust standard deviation:',result[7])
        print('>>> 0.1586553 percentile.....:',result[8])
        print('>>> 0.8413447 percentile.....:',result[9])

    return result

#------------------------------------------------------------------------------


