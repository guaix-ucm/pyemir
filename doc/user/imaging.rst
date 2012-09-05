
Imaging Recipes
===============

Stare Image
-----------

The effect of recording images of the sky in a given pointing 
position of the TS.

Procedure
+++++++++

Nodded/Beam-switched images
---------------------------

The effect of recording a series of stare images, with the same
acquisition parameters, and taken by pointing the TS in cycles
between two, or more, sky positions. Displacements are larger
than the EMIR FOV, so the images have no common area. Used
for sky subtraction.

Procedure
+++++++++

Dithered images
---------------
The effect of recording a series of stare images, with the same
acquisition parameters, and taken by pointing to a number of
sky positions, with separations of the order of arcsec, either by
nodding the TS, tilting the TS M2 or shifting the EMIR DTU.
Displacements are of the order of several pixels (even
fractional). Images share the large majority of the sky positions
so they can be coadded. Used for avoid cosmetic effects and/or
improve the SNR. Superflat and/or supersky frames can be built
from the image series.

Procedure
+++++++++
Images are corrected from dark, non-linearity and flat. Then, an iterative
process starts

Micro-dithered images
---------------------
The effect of recording a series of stare images, with the same
acquisition parameters, and taken by pointing to a number of
sky positions, with separations of the order of sub arcsecs,
either by moving the either by nodding the TS, tilting the TS
M2 or shifting the EMIR DTU, the latter being the most likely
option. Displacements are of the order of fraction of pixels.
Images share the large majority of the sky positions so they can
be coadded. Used for improving the spatial resolution of the
resulting images and not valid for sky or superflat images.

Procedure
+++++++++

