
Auxiliary Recipes
=================

Bias Image
----------
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

Products
++++++++

+-------------------+-------------------------------------------------------+
| Name              | Type                                                  |
+===================+=======================================================+
| ``'darkframe'``   | :class:`~emir.dataproducts.MasterDark`                | 
+-------------------+-------------------------------------------------------+
| ``'stats'``       | :class:`~emir.dataproducts.ChannelLevelStatistics`    |
+-------------------+-------------------------------------------------------+

Procedure
+++++++++

Intensity Flat-Field
--------------------
The required actions to set the TS and EMIR at the
configuration from which sky and/or artificial illumination flat
field data acquisition can proceed and take data.

Procedure
+++++++++

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

Products
++++++++

+-------------------+-------------------------------------------------------+
| Name              | Type                                                  |
+===================+=======================================================+
| ``'flatframe'``   | :class:`~emir.dataproducts.MasterIntensityFlat`       | 
+-------------------+-------------------------------------------------------+



