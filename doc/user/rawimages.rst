
Raw Image data products
=======================

Readout modes
*************

Single
------
In single mode, the detector is readout after a reset. After the readout, the
detector is reset gain. This readout mode has limited utillity and it is meant
for engeneering.

The frame contains a single HDU. The data section of the HDU contains
the 2048x2048 dataframe. In the header, the following keywords are set

 ::

   READPROC = F / The frame has been preprocessed after readout
   READMODE = ‘SINGLE’  / The readmode used to adquire the image


Correlated double sampling
--------------------------
In correlated double sampling, the detector is reset, then read inmediately
after the reset and then read after the programed exposure time.

The frame contains a single HDU. The data section of the HDU contains
a 2048x2048x2 dataframe. The first layer contains the readout after reset
and the second the readout after the exposure time.  In the header, the 
following keywords are set

 ::

   READPROC = F / The frame has been preprocessed after readout
   READMODE = ‘CDS’  / The readmode used to adquire the image

Fowler
------
Fowler mode is an extension of CDS. The detector is reset, then read ``n``
times, exposed and then read agin ``n`` times. The exposure time in this
case is equal to the time between correlated reads.

The frame contains a single HDU. The data section of the HDU contains
a 2048x2048x2n dataframe. The first ``n`` layers contain the readouts after 
reset and the second ``n`` the readouts after the exposure time.  
In the header, the following keywords are set

 ::

   READPROC = F / The frame has been preprocessed after readout
   READMODE = ‘FOWLER’  / The readmode used to adquire the image
   READSAMP = 2*n / Number of samples taken

Follow-up-the ramp
------------------
In ramp mode, the exposure time is sampled ``n`` times.

The frame contains a single HDU. The data section of the HDU contains 
a 2048x2048xn dataframe. Each layer contains the n-th readouts after 
reset.
In the header, the following keywords are set

 ::

    READPROC = F / The frame has been preprocessed after readout
    READMODE = ‘RAMP’  / The readmode used to adquire the image
    READSAMP = n / Number of samples taken


Image types
***********

Bias
----

Dark
----

Flat
----

Target
------

