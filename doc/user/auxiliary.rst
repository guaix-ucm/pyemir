
Auxiliary Recipes
=================

Bias Image
----------

:Mode: Bias Image
:Recipe class: :class:`~emir.recipes.BiasRecipe`
:Input class: :class:`~emir.recipes.BiasRecipeInput`
:Result class: :class:`~emir.recipes.BiasRecipeResult`

The actions to calibrate the zero (pedestal) level of the detector
plus associated control electronic by taken images with null
integration time.

Requeriments
++++++++++++

+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'master_bpm'``         | Product       | NA         |      Master BPM frame         |
+--------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++
The frames in the observed block are stacked together using the median of them as the final result.
The variance of the result frame is computed using two different methods. 
The first method computes the variance across the pixels in the different frames stacked.
The second method computes the variance en each channel in the result frame.

Products
++++++++

+-------------------+-------------------------------------------------------+
| Name              | Type                                                  |
+===================+=======================================================+
| ``'biasframe'``   | :class:`~emir.dataproducts.MasterBias`                | 
+-------------------+-------------------------------------------------------+
| ``'stats'``       | :class:`~emir.dataproducts.ChannelLevelStatistics`    |
+-------------------+-------------------------------------------------------+

Dark Current Image
------------------

:Mode: Dark Current Image
:Recipe class: :class:`~emir.recipes.DarkRecipe`
:Input class: :class:`~emir.recipes.DarkRecipeInput`
:Result class: :class:`~emir.recipes.DarkRecipeResult`

The actions to measure the variation of the intrinsic signal of the
system by taken images under zero illumination condition and
long integration time.

Requeriments
++++++++++++

+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'master_bpm'``         | Product       | NA         |      Master BPM frame         |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``        | Product       | NA         | Master Bias frame             |
+--------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++
The frames in the observed block are subtracted from the master bias and then they are stacked together, using the median of them as the final result.
The variance of the result frame is computed using two different methods. 
The first method computes the variance across the pixels in the different frames stacked.
The second method computes the variance en each channel in the result frame.

Products
++++++++

+-------------------+-------------------------------------------------------+
| Name              | Type                                                  |
+===================+=======================================================+
| ``'darkframe'``   | :class:`~emir.dataproducts.MasterDark`                | 
+-------------------+-------------------------------------------------------+
| ``'stats'``       | :class:`~emir.dataproducts.ChannelLevelStatistics`    |
+-------------------+-------------------------------------------------------+

.. _ff-recipe-label:

Intensity Flat-Field
--------------------

:Mode: Intensity Flat-Field
:Recipe class: :class:`~emir.recipes.IntensityFlatRecipe`
:Input class: :class:`~emir.recipes.IntensityFlatRecipeInput`
:Result class: :class:`~emir.recipes.IntensityFlatRecipeResult`

The required actions to set the TS and EMIR at the
configuration from which sky and/or artificial illumination flat
field data acquisition can proceed and take data.

Requeriments
++++++++++++

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

Procedure
+++++++++
The frames in the observed block are subtracted from the master bias and the master dark.
The frames are corrected from non-linearity.

The frames with lamps-on and with lamps-off are stacked using the median, and then the
combined lamps-off frame is subtracted from the lamps-on frame. The result is the
subtracted frame, scaled to have a mean value of 1.


Products
++++++++

+-------------------+-------------------------------------------------------+
| Name              | Type                                                  |
+===================+=======================================================+
| ``'flatframe'``   | :class:`~emir.dataproducts.MasterIntensityFlat`       | 
+-------------------+-------------------------------------------------------+



