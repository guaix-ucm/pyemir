
Imaging Recipes
===============

Stare Image
-----------

The effect of recording images of the sky in a given pointing 
position of the TS.

Inputs
++++++

+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'master_bpm'``         | Product       | NA         |      Master BPM frame         |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``        | Product       | NA         | Master Bias frame             |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_dark'``        | Product       | NA         | Master Dark frame             |
+--------------------------+---------------+------------+-------------------------------+
| ``'nonlinearity'``       | Product       | [1.0, 0.0] | Master non-linearity          |
|                          |               |            | calibration                   |
+--------------------------+---------------+------------+-------------------------------+
|``'master_intensity_ff'`` | Product       | NA         | Master Intensity flat-field   |
|                          |               |            | frame                         |
+--------------------------+---------------+------------+-------------------------------+
| ``'extinction '``        | Parameter     | 0.0        | Mean atmospheric extinction   |
+--------------------------+---------------+------------+-------------------------------+
| ``'sources    '``        | Parameter     | None       | List of (x, y) coordinates to |
|                          |               |            | measure FWHM                  |
+--------------------------+---------------+------------+-------------------------------+
| ``'offsets    '``        | Parameter     | None       | List of pairs of offsets      |
+--------------------------+---------------+------------+-------------------------------+
| ``'iterations'``         | Parameter     | 4          | Iterations of the recipe      |
+--------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++


Results
+++++++

The result of the Recipe is an object of type :class:`emir.recipes.StareImageRecipeResult`. 
It contains two objects, a :class:`emir.dataproducts.FrameDataProduct` containing the result frame
and a :class:`emir.dataproducts.SourcesCatalog` containing a catalog of sources.

Nodded/Beam-switched images
---------------------------

The effect of recording a series of stare images, with the same
acquisition parameters, and taken by pointing the TS in cycles
between two, or more, sky positions. Displacements are larger
than the EMIR FOV, so the images have no common area. Used
for sky subtraction.

Inputs
++++++

+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'master_bpm'``         | Product       | NA         |      Master BPM frame         |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``        | Product       | NA         | Master Bias frame             |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_dark'``        | Product       | NA         | Master Dark frame             |
+--------------------------+---------------+------------+-------------------------------+
| ``'nonlinearity'``       | Product       | [1.0, 0.0] | Master non-linearity          |
|                          |               |            | calibration                   |
+--------------------------+---------------+------------+-------------------------------+
|``'master_intensity_ff'`` | Product       | NA         | Master Intensity flat-field   |
|                          |               |            | frame                         |
+--------------------------+---------------+------------+-------------------------------+
| ``'extinction '``        | Parameter     | 0.0        | Mean atmospheric extinction   |
+--------------------------+---------------+------------+-------------------------------+
| ``'sources    '``        | Parameter     | None       | List of (x, y) coordinates to |
|                          |               |            | measure FWHM                  |
+--------------------------+---------------+------------+-------------------------------+
| ``'offsets    '``        | Parameter     | None       | List of pairs of offsets      |
+--------------------------+---------------+------------+-------------------------------+
| ``'iterations'``         | Parameter     | 4          | Iterations of the recipe      |
+--------------------------+---------------+------------+-------------------------------+



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

Inputs
++++++

+------------------------------+---------------+------------------+-------------------------------+
| Name                         | Type          | Default          | Meaning                       |
+==============================+===============+==================+===============================+
| ``'master_bpm'``             | Product       | NA               |      Master BPM frame         |
+------------------------------+---------------+------------------+-------------------------------+
| ``'master_bias'``            | Product       | NA               | Master Bias frame             |
+------------------------------+---------------+------------------+-------------------------------+
| ``'master_dark'``            | Product       | NA               | Master Dark frame             |
+------------------------------+---------------+------------------+-------------------------------+
| ``'nonlinearity'``           | Product       | [1.0, 0.0]       | Master non-linearity          |
|                              |               |                  | calibration                   |
+------------------------------+---------------+------------------+-------------------------------+
|``'master_intensity_ff'``     | Product       | NA               | Master Intensity flat-field   |
|                              |               |                  | frame                         |
+------------------------------+---------------+------------------+-------------------------------+
| ``'extinction'``             | Parameter     | 0.0              | Mean atmospheric extinction   |
+------------------------------+---------------+------------------+-------------------------------+
| ``'sources'``                | Parameter     | None             | List of (x, y) coordinates to |
|                              |               |                  | measure FWHM                  |
+------------------------------+---------------+------------------+-------------------------------+
| ``'offsets'``                | Parameter     | None             | List of pairs of offsets      |
+------------------------------+---------------+------------------+-------------------------------+
| ``'iterations'``             | Parameter     | 4                | Iterations of the recipe      |
+------------------------------+---------------+------------------+-------------------------------+
| ``'sky_images'``             | Parameter     | 5                | Images used to estimate the   | 
|                              |               |                  | background before and after   |
|                              |               |                  | current image                 |
+------------------------------+---------------+------------------+-------------------------------+
| ``'sky_images_sep_time'``    | Parameter     | 10               | Maximum separation time       |
|                              |               |                  | between consecutive sky images| 
|                              |               |                  | in minutes                    |
+------------------------------+---------------+------------------+-------------------------------+
|``'check_photometry_levels'`` | Parameter     | [0.5, 0.8]       | Levels to check the flux of   |
|                              |               |                  | the objects                   |
+------------------------------+---------------+------------------+-------------------------------+
|``'chec_photometry_actions'`` | Parameter     | ['warn', 'warn', | Actions to take on images     |
|                              |               | 'default']       |                               |     
+------------------------------+---------------+------------------+-------------------------------+


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

Inputs
++++++

+------------------------------+---------------+------------------+-------------------------------+
| Name                         | Type          | Default          | Meaning                       |
+==============================+===============+==================+===============================+
| ``'master_bpm'``             | Product       | NA               |      Master BPM frame         |
+------------------------------+---------------+------------------+-------------------------------+
| ``'master_bias'``            | Product       | NA               | Master Bias frame             |
+------------------------------+---------------+------------------+-------------------------------+
| ``'master_dark'``            | Product       | NA               | Master Dark frame             |
+------------------------------+---------------+------------------+-------------------------------+
| ``'nonlinearity'``           | Product       | [1.0, 0.0]       | Master non-linearity          |
|                              |               |                  | calibration                   |
+------------------------------+---------------+------------------+-------------------------------+
|``'master_intensity_ff'``     | Product       | NA               | Master Intensity flat-field   |
|                              |               |                  | frame                         |
+------------------------------+---------------+------------------+-------------------------------+
| ``'extinction'``             | Parameter     | 0.0              | Mean atmospheric extinction   |
+------------------------------+---------------+------------------+-------------------------------+
| ``'sources'``                | Parameter     | None             | List of (x, y) coordinates to |
|                              |               |                  | measure FWHM                  |
+------------------------------+---------------+------------------+-------------------------------+
| ``'offsets'``                | Parameter     | None             | List of pairs of offsets      |
+------------------------------+---------------+------------------+-------------------------------+
| ``'iterations'``             | Parameter     | 4                | Iterations of the recipe      |
+------------------------------+---------------+------------------+-------------------------------+
| ``'sky_images'``             | Parameter     | 5                | Images used to estimate the   | 
|                              |               |                  | background before and after   |
|                              |               |                  | current image                 |
+------------------------------+---------------+------------------+-------------------------------+
| ``'sky_images_sep_time'``    | Parameter     | 10               | Maximum separation time       |
|                              |               |                  | between consecutive sky images| 
|                              |               |                  | in minutes                    |
+------------------------------+---------------+------------------+-------------------------------+
|``'check_photometry_levels'`` | Parameter     | [0.5, 0.8]       | Levels to check the flux of   |
|                              |               |                  | the objects                   |
+------------------------------+---------------+------------------+-------------------------------+
|``'chec_photometry_actions'`` | Parameter     | ['warn', 'warn', | Actions to take on images     |
|                              |               | 'default']       |                               |     
+------------------------------+---------------+------------------+-------------------------------+
|``'subpixelization'``         | Parameter     | 4                | Number of subdivision of each |
|                              |               |                  | pixel side                    |     
+------------------------------+---------------+------------------+-------------------------------+
|``'window'``                  | Parameter     | None             | Region of interest            |
+------------------------------+---------------+------------------+-------------------------------+


Procedure
+++++++++

