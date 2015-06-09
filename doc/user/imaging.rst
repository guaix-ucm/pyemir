
Imaging Recipes
===============

Stare Image
-----------

:Mode: Stare Image
:Recipe class: :class:`~emirdrp.recipes.StareImageRecipe`
:Input class: :class:`~emirdrp.recipes.StareImageRecipeInput`
:Result class: :class:`~emirdrp.recipes.StareImageRecipeResult`

The effect of recording images of the sky in a given pointing 
position of the TS.

Requeriments
++++++++++++

=========================== ========== =========== ==============================
 Name                       Type       Default     Meaning                       
=========================== ========== =========== ==============================
``'master_bpm'``            Product     NA         Master BPM frame              
``'master_bias'``           Product     NA         Master Bias frame             
``'master_dark'``           Product     NA         Master Dark frame             
``'nonlinearity'``          Product     [1.0, 0.0] Master non-linearity          
                                                   calibration                   
``'master_intensity_ff'``   Product     NA         Master Intensity flat-field   
                                                   frame                         
``'extinction'``            Parameter   0.0        Mean atmospheric extinction   
``'sources'``               Parameter   None       List of (x, y) coordinates to 
                                                   measure FWHM                  
``'offsets'``               Parameter   None       List of pairs of offsets      
``'iterations'``            Parameter   4          Iterations of the recipe      
=========================== ========== =========== ==============================

Procedure
+++++++++

The block of raw frames are processed in several stages. First, the frames
are corrected from bias, dark, non-linearity and intensity flat field. For
these steps we use the frames in the requeriments from ``'master_bpm'``
to ``'master_intensity_ff'``

.. note::
   It is not clear the need to have a master_bias requirement, as it is
   only needed in one readout mode

If the parameter ``offsets`` is ``None``, the recipe computes
offset information using the WCS stored in each frame FITS header. If it
is defined, the information in ``offsets`` is used instead.

The size of the result frame is computed using the sizes of the input
frames and the offstes between them. New intermediate frames are
created resizing the input frames.

A sky flat is computed using the input frames. Each frame is scaled
acording to its mean value and then they are combined together
using a sigma clipping algorithm with low and high rejection limits
of 3 sigma.  Each input frame is then divided by the sky flat. 

Next, the sky level is estimated in each frame by obtaining the median.
The sky value is then subtracted of each frame. The frames 
(sky-subtracted, flat-fielded and resized) are then stacked. The frames
are scaled acording to its airmass and the value of ``'extinction'``.
The algorithm used to stack frames is :func:`numina.array.combine.quantileclip`
with 10% points rejected at both ends of the distribution.

Results
+++++++
The result of the Recipe is an object of type :class:`~emirdrp.recipes.StareImageRecipeResult`. 
It contains two objects, a :class:`~emirdrp.dataproducts.FrameDataProduct` containing the result frame
and a :class:`~emirdrp.dataproducts.SourcesCatalog` containing a catalog of sources.

Nodded/Beam-switched images
---------------------------

:Mode: Nodded/Beam-switched images
:Recipe class: :class:`~emirdrp.recipes.NBImageRecipe`
:Input class: :class:`~emirdrp.recipes.NBImageRecipeInput`
:Result class: :class:`~emirdrp.recipes.NBImageRecipeResult`

The effect of recording a series of stare images, with the same
acquisition parameters, and taken by pointing the TS in cycles
between two, or more, sky positions. Displacements are larger
than the EMIR FOV, so the images have no common area. Used
for sky subtraction.

Requeriments
++++++++++++

+---------------------------+---------------+------------+-------------------------------+
| Name                      | Type          | Default    | Meaning                       |
+===========================+===============+============+===============================+
| ``'master_bpm'``          | Product       | NA         |      Master BPM frame         |
+---------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``         | Product       | NA         | Master Bias frame             |
+---------------------------+---------------+------------+-------------------------------+
| ``'master_dark'``         | Product       | NA         | Master Dark frame             |
+---------------------------+---------------+------------+-------------------------------+
| ``'nonlinearity'``        | Product       | [1.0, 0.0] | Master non-linearity          |
|                           |               |            | calibration                   |
+---------------------------+---------------+------------+-------------------------------+
| ``'master_intensity_ff'`` | Product       | NA         | Master Intensity flat-field   |
|                           |               |            | frame                         |
+---------------------------+---------------+------------+-------------------------------+
| ``'extinction'``          | Parameter     | 0.0        | Mean atmospheric extinction   |
+---------------------------+---------------+------------+-------------------------------+
| ``'sources'``             | Parameter     | None       | List of (x, y) coordinates to |
|                           |               |            | measure FWHM                  |
+---------------------------+---------------+------------+-------------------------------+
| ``'offsets'``             | Parameter     | None       | List of pairs of offsets      |
+---------------------------+---------------+------------+-------------------------------+
| ``'iterations'``          | Parameter     | 4          | Iterations of the recipe      |
+---------------------------+---------------+------------+-------------------------------+



Procedure
+++++++++
The block of raw frames contains both sky and target images. They are treated differently at some
stages. Sky frames have ``IMGTYP = 'SKY'`` in their FITS headers. Target frames have 
``IMGTYP = 'TARGET'``. 

All the frames are corrected from bias, dark, non-linearity and intensity flat field. For
these steps we use the frames in the requeriments from ``'master_bpm'``
to ``'master_intensity_ff'``

.. note::
   It is not clear the need to have a master_bias requirement, as it is
   only needed in one readout mode

Then, an iterative process starts. The number of iterations is controlled by the
parameter ``'iterations'``.

Base step
'''''''''
Offsets between the target frames are obtained. If the parameter ``offsets`` 
is ``None``, the recipe computes
offset information using the WCS stored in each frame FITS header. If it
is defined, the information in ``offsets`` is used instead.

The size of the result frame is computed using the sizes of the target
frames and the offsets between them. New intermediate frames are
created resizing the target input frames.

A sky flat is computed using the input sky frames. Each sky frame is scaled
acording to its mean value and then they are combined together
using a sigma clipping algorithm with low and high rejection limits
of 3 sigma.  Each input target frame is then divided by the sky flat. 

Next, the sky level is estimated in each frame by obtaining the median of the
nearest sky image.
The sky value is then subtracted of each frame. The target frames 
(sky-subtracted, flat-fielded and resized) are then stacked. The frames
are scaled acording to its airmass and the value of ``'extinction'``.
The algorithm used to stack frames is :func:`numina.array.combine.quantileclip`
with 10% points rejected at both ends of the distribution.

Check step
''''''''''
In the next step, several checkings are performed in the result image.

The centroids of bright objects are compared between the input target
frames and the result frame. This test allows to check if the
offsets are correct and to refine the offsets.

The flux of bright objects is compared between the input target frames
and the result frame. This test allows to find target frames with
abnormal illumination (due to clouds, for example). The 
parameter ``'check_photometry_levels'`` mark different categories
of clasification of the frames acording the fraction of the median
flux level of the frames. The parameter ``'check_photometry_actions'``
allow the user to select the action to take in each category.
The allowed actions are ``'default'`, ``'warn'`` and ``'reject'``.

.. warning::
   The offset-recompute routine is not yet implemented

Full reduction step
'''''''''''''''''''
Using the latest available result image (in the first iteration, that of the base step), 
a segmentation mask is computed. This segmentation mask applies to target frames only.

.. note::
   A segmentation mask for each **sky frame** is being considered

The sky flat is applied to the target frames.

The sky level for target frames is estimated using the median value of the nearest
sky frames in the observed series. We use a number of 
``'sky_images'`` frames before and after and never separated more than 
``'sky_images_sep_time'`` minutes.

The target frames (sky-subtracted, flat-fielded and resized) are then stacked. The frames
are scaled acording to its airmass and the value of ``'extinction'``.
The algorithm used to stack frames is :func:`numina.array.combine.quantileclip`
with 10% points rejected at both ends of the distribution.

This last step is repeated ``'iterations'`` times, the segmentation mask computed
from the result of the previous step.

Results
+++++++
The result of the Recipe is an object of type :class:`~emirdrp.recipes.NBImageRecipeResult`. 
It contains two objects, a :class:`~emirdrp.dataproducts.FrameDataProduct` containing the result frame
and a :class:`~emirdrp.dataproducts.SourcesCatalog` containing a catalog of sources.

Dithered images
---------------

:Mode: Dithered images
:Recipe class: :class:`~emirdrp.recipes.DitheredImageRecipe`
:Input class: :class:`~emirdrp.recipes.DitheredImageRecipeInput`
:Result class: :class:`~emirdrp.recipes.DitheredImageRecipeResult`

The effect of recording a series of stare images, with the same
acquisition parameters, and taken by pointing to a number of
sky positions, with separations of the order of arcsec, either by
nodding the TS, tilting the TS M2 or shifting the EMIR DTU.
Displacements are of the order of several pixels (even
fractional). Images share the large majority of the sky positions
so they can be coadded. Used for avoid cosmetic effects and/or
improve the SNR. Superflat and/or supersky frames can be built
from the image series.

Requeriments
++++++++++++
+-------------------------------+---------------+------------------+-------------------------------+
| Name                          | Type          | Default          | Meaning                       |
+===============================+===============+==================+===============================+
| ``'master_bpm'``              | Product       | NA               |      Master BPM frame         |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'master_bias'``             | Product       | NA               | Master Bias frame             |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'master_dark'``             | Product       | NA               | Master Dark frame             |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'nonlinearity'``            | Product       | [1.0, 0.0]       | Master non-linearity          |
|                               |               |                  | calibration                   |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'master_intensity_ff'``     | Product       | NA               | Master Intensity flat-field   |
|                               |               |                  | frame                         |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'extinction'``              | Parameter     | 0.0              | Mean atmospheric extinction   |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'sources'``                 | Parameter     | None             | List of (x, y) coordinates to |
|                               |               |                  | measure FWHM                  |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'offsets'``                 | Parameter     | None             | List of pairs of offsets      |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'iterations'``              | Parameter     | 4                | Iterations of the recipe      |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'sky_images'``              | Parameter     | 5                | Images used to estimate the   | 
|                               |               |                  | background before and after   |
|                               |               |                  | current image                 |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'sky_images_sep_time'``     | Parameter     | 10               | Maximum separation time       |
|                               |               |                  | between consecutive sky images| 
|                               |               |                  | in minutes                    |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'check_photometry_levels'`` | Parameter     | [0.5, 0.8]       | Levels to check the flux of   |
|                               |               |                  | the objects                   |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'chec_photometry_actions'`` | Parameter     | ['warn', 'warn', | Actions to take on images     |
|                               |               | 'default']       |                               |
+-------------------------------+---------------+------------------+-------------------------------+


Procedure
+++++++++

The block of raw frames are processed in several stages. First, the frames
are corrected from bias, dark, non-linearity and intensity flat field. For
these steps we use the frames in the requeriments from ``'master_bpm'``
to ``'master_intensity_ff'``

.. note::
   It is not clear the need to have a master_bias requirement, as it is
   only needed in one readout mode

Then, an iterative process starts. The number of iterations is controlled by the
parameter ``'iterations'``.

Base step
'''''''''
Offsets between the frames are obtained. If the parameter ``offsets`` 
is ``None``, the recipe computes
offset information using the WCS stored in each frame FITS header. If it
is defined, the information in ``offsets`` is used instead.

The size of the result frame is computed using the sizes of the input
frames and the offstes between them. New intermediate frames are
created resizing the input frames.

A sky flat is computed using the input frames. Each frame is scaled
acording to its mean value and then they are combined together
using a sigma clipping algorithm with low and high rejection limits
of 3 sigma.  Each input frame is then divided by the sky flat. 

Next, the sky level is estimated in each frame by obtaining the median.
The sky value is then subtracted of each frame. The frames 
(sky-subtracted, flat-fielded and resized) are then stacked. The frames
are scaled acording to its airmass and the value of ``'extinction'``.
The algorithm used to stack frames is :func:`numina.array.combine.quantileclip`
with 10% points rejected at both ends of the distribution.

Check step
''''''''''
In the next step, several checkings are performed in the result image.

The centroids of bright objects are compared between the input
frames and the result frame. This test allows to check if the
offsets are correct and to refine the offsets.

The flux of bright objects is compared between the input frames
and the result frame. This test allows to find frames with
abnormal illumination (due to clouds, for eample). The 
parameter ``'check_photometry_levels'`` mark different categories
of clasification of the frames acording the fraction of the median
flux level of the frames. The parameter ``'check_photometry_actions'``
allow the user to select the action to take in each category.
The allowed actions are ``'default'`, ``'warn'`` and ``'reject'``.

.. warning::
   The offset-recompute routine is not yet implemented

Full reduction step
'''''''''''''''''''
Using the latest available result image (in the first iteration, that of the base step), 
a segmentation mask is computed.
The segmentation mask is used to avoid objects when computing a new sky flat.
With the frames corrected with the new sky flat, the sky level is estimated.
For each frame, we use frames before and after in the series to compute a
median sky, that is subtracted from each frame. We use a number of 
``'sky_images'`` frames before and after and never separated more than 
``'sky_images_sep_time'`` minutes.

The frames (sky-subtracted, flat-fielded and resized) are then stacked. The frames
are scaled acording to its airmass and the value of ``'extinction'``.
The algorithm used to stack frames is :func:`numina.array.combine.quantileclip`
with 10% points rejected at both ends of the distribution.

This last step is repeated ``'iterations'`` times, the segmentation mask computed
from the result of the previous step.


Results
+++++++
The result of the Recipe is an object of type :class:`~emirdrp.recipes.DitheredImageRecipeResult`. 
It contains two objects, a :class:`~emirdrp.dataproducts.FrameDataProduct` containing the result frame
and a :class:`~emirdrp.dataproducts.SourcesCatalog` containing a catalog of sources.

Micro-dithered images
---------------------

:Mode: Micro-dithered images
:Recipe class: :class:`~emirdrp.recipes.MicroDitheredImageRecipe`
:Input class: :class:`~emirdrp.recipes.MicroDitheredImageRecipeInput`
:Result class: :class:`~emirdrp.recipes.MicroDitheredImageRecipeResult`

The effect of recording a series of stare images, with the same
acquisition parameters, and taken by pointing to a number of
sky positions, with separations of the order of sub arcsecs,
either by moving the either by nodding the TS, tilting the TS
M2 or shifting the EMIR DTU, the latter being the most likely
option. Displacements are of the order of fraction of pixels.
Images share the large majority of the sky positions so they can
be coadded. Used for improving the spatial resolution of the
resulting images and not valid for sky or superflat images.


Requeriments
++++++++++++

+-------------------------------+---------------+------------------+-------------------------------+
| Name                          | Type          | Default          | Meaning                       |
+===============================+===============+==================+===============================+
| ``'master_bpm'``              | Product       | NA               |      Master BPM frame         |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'master_bias'``             | Product       | NA               | Master Bias frame             |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'master_dark'``             | Product       | NA               | Master Dark frame             |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'nonlinearity'``            | Product       | [1.0, 0.0]       | Master non-linearity          |
|                               |               |                  | calibration                   |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'master_intensity_ff'``     | Product       | NA               | Master Intensity flat-field   |
|                               |               |                  | frame                         |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'extinction'``              | Parameter     | 0.0              | Mean atmospheric extinction   |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'sources'``                 | Parameter     | None             | List of (x, y) coordinates to |
|                               |               |                  | measure FWHM                  |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'offsets'``                 | Parameter     | None             | List of pairs of offsets      |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'iterations'``              | Parameter     | 4                | Iterations of the recipe      |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'sky_images'``              | Parameter     | 5                | Images used to estimate the   | 
|                               |               |                  | background before and after   |
|                               |               |                  | current image                 |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'sky_images_sep_time'``     | Parameter     | 10               | Maximum separation time       |
|                               |               |                  | between consecutive sky images| 
|                               |               |                  | in minutes                    |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'check_photometry_levels'`` | Parameter     | [0.5, 0.8]       | Levels to check the flux of   |
|                               |               |                  | the objects                   |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'chec_photometry_actions'`` | Parameter     | ['warn', 'warn', | Actions to take on images     |
|                               |               | 'default']       |                               |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'subpixelization'``         | Parameter     | 4                | Number of subdivision of each |
|                               |               |                  | pixel side                    |
+-------------------------------+---------------+------------------+-------------------------------+
| ``'window'``                  | Parameter     | None             | Region of interest            |
+-------------------------------+---------------+------------------+-------------------------------+


Procedure
+++++++++

The procedure followed by this recipe is equivalent to Dithered images. They differ in the aspects
controlled by the parameters ``'subpixelization'`` and  ``'window'``. If ``'window'`` is different
to ``None``, the frames are clipped to the size ``'window'``. Each pixel of the input frames
is subdivided in ``'subpixelization'`` x ``'subpixelization'`` pixels. 

Results
+++++++
The result of the Recipe is an object of type :class:`~emirdrp.recipes.MicroDitheredImageRecipeResult`. 
It contains two objects, a :class:`~emirdrp.dataproducts.FrameDataProduct` containing the result frame
and a :class:`~emirdrp.dataproducts.SourcesCatalog` containing a catalog of sources.

