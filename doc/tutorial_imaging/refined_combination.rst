===============================
Improving the image combination
===============================

The image combination can be improve by tuning some of the parameters of the
recipe ``FULL_DITHERED_IMAGE`` (step 2 in the previous section).
In this sense, there is no need to repeat the basic reduction of the individual
exposures (step1).

As previously mentioned, two are the problems that we want to solve:

1. **Improve the offsets between individual exposures:** this can be achieved
   in several ways:

   - by setting the requirement ``refine_offsets: True``: in this case a
     cross-correlation between subimage regions around bright targets is
     carried out to derive refined offsets. See subsection
     :ref:`improving_offsets_1` below.

   - by providing an ASCII file with a list of offsets measured independently 
     by the user and indicated with the requirement ``offsets:
     user_offsets.txt``. See subsection :ref:`improving_offsets_2` below.

   - by providing the same ASCII file with precomputed offsets (as in the
     previous item) and using, in addition, the cross-correlation method. In
     this case, both requirements ``refine_offsets: True`` and ``offsets:
     user_offsets.txt`` must be set. See subsection :ref:`improving_offsets_3`
     below.

2. **Improve the sky background level estimation:** the background level can be
   improved by:

   - generating an object mask and iterating the combination process. See
     subsection :ref:`improving_skybackground_1` below.

   - introducing an *ad hoc* fit to a low-order polynomial surface to the sky
     background. See subsection :ref:`improving_skybackground_2` below.


.. _improving_offsets_1:

Improving offsets (method #1)
-----------------------------

We can activate the use of 2D cross-correlation of subimages around bright
targets to obtain refined offsets. This method works only if the initial
offsets (either derived from the WCS information in the image headers or from
an external file provided by the user) are a reasonable approximation to the
refined values. To activate this option it it necessary to set the requirement
``refine_offsets: True`` in the observation result file.

This option is already set in line number 120 of the file ``dithered_v1.yaml``,
available in the downloaded package for this tutorial.

.. literalinclude:: dithered_v1.yaml
   :lines: 117-121
   :emphasize-lines: 4
   :linenos:
   :lineno-start: 117

The refined version of the combined image is then obtained by executing numina
again with this new observation result file:

::

   (emir) $ numina run dithered_v1.yaml --link-files -r control.yaml

.. generada con --geometry 0,0,850,1200
.. convert combined_v1.png -trim combined_v1_trimmed.png
.. convert -delay 100 -loop 0 combined_v[01]_trimmed.png comparison_v1.gif

.. only:: html

   .. image:: comparison_v1.gif
      :width: 100%
      :alt: combined image, version 1 compared with version 0

.. only:: latex

   |pic_combined_v0_trimmed| |pic_combined_v1_trimmed|

   .. |pic_combined_v0_trimmed| image:: combined_v0_trimmed.png
      :width: 48%

   .. |pic_combined_v1_trimmed| image:: combined_v1_trimmed.png
      :width: 48%

.. generada con --geometry 0,0,850,1200 --bbox 1100,1600,800,1300
.. convert combined_v1_zoom.png -trim combined_v1_zoom_trimmed.png
.. convert -delay 100 -loop 0 combined_v[01]_zoom_trimmed.png comparison_v1_zoom.gif

.. only:: html

   .. image:: comparison_v1_zoom.gif
      :width: 100%
      :alt: combined image, version 1 compared with version 0

.. only:: latex

   |pic_combined_v0_zoom_trimmed| |pic_combined_v1_zoom_trimmed|

   .. |pic_combined_v0_zoom_trimmed| image:: combined_v0_zoom_trimmed.png
      :width: 48%

   .. |pic_combined_v1_zoom_trimmed| image:: combined_v1_zoom_trimmed.png
      :width: 48%

.. _improving_offsets_2:

Improving offsets (method #2)
-----------------------------

An alternative to the use of the offsets computed from the WCS information in
the image header is to provide a two-column ASCII file with the measured
offsets between the individual images. The (arbitray) name of that file must be
provided through the requirement ``offsets:``. For this tutoral, we are
providing such a file with the name ``user_offsets.txt``. Note that this file
must be placed within the ``data`` subdirectory.

The observation result file ``dithered_v2.yaml`` is similar to the initial
``dithered_v0.yaml`` file, with the inclusion of the new requirement (line
number 121):

.. literalinclude:: dithered_v2.yaml
   :lines: 117-122
   :emphasize-lines: 4
   :linenos:
   :lineno-start: 117

The contents of the ASCII file with the measured offsets is the following:

::

   (emir) cat data/user_offsets.txt
   822 907
   730 660
   555 863
   620 998
   895 741
   545 674
   708 811
   830 911
   735 666
   561 868
   626 1003
   901 746
   551 679
   715 816

Execute numina to obtain the new version of the combined image:

::

   (emir) $ numina run dithered_v2.yaml --link-files -r control.yaml

.. generada con --geometry 0,0,850,1200
.. convert combined_v2.png -trim combined_v2_trimmed.png
.. convert -delay 100 -loop 0 combined_v[02]_trimmed.png comparison_v2.gif

.. only:: html

   .. image:: comparison_v2.gif
      :width: 100%
      :alt: combined image, version 2 compared with version 0

.. only:: latex

   |pic_combined_v0_trimmed| |pic_combined_v2_trimmed|

   .. |pic_combined_v0_trimmed| image:: combined_v0_trimmed.png
      :width: 48%

   .. |pic_combined_v2_trimmed| image:: combined_v2_trimmed.png
      :width: 48%

.. generada con --geometry 0,0,850,1200 --bbox 1100,1600,800,1300
.. convert combined_v2_zoom.png -trim combined_v2_zoom_trimmed.png
.. convert -delay 100 -loop 0 combined_v[02]_zoom_trimmed.png comparison_v2_zoom.gif

.. only:: html

   .. image:: comparison_v2_zoom.gif
      :width: 100%
      :alt: combined image, version 2 compared with version 0

.. only:: latex

   |pic_combined_v0_zoom_trimmed| |pic_combined_v2_zoom_trimmed|

   .. |pic_combined_v0_zoom_trimmed| image:: combined_v0_zoom_trimmed.png
      :width: 48%

   .. |pic_combined_v2_zoom_trimmed| image:: combined_v2_zoom_trimmed.png
      :width: 48%



.. _improving_offsets_3:

Improving offsets (method #3)
-----------------------------

It is also possible to combine both the offsets provided by the user through an
external ASCII file, as well as the cross-correlation method to improve those
numbers.

The last lines of the new observation result file ``dithered_v3.yaml`` are the
following:

.. literalinclude:: dithered_v3.yaml
   :lines: 117-122
   :emphasize-lines: 4-5
   :linenos:
   :lineno-start: 117

Execute numina again with this new observation result file:

::

   (emir) $ numina run dithered_v3.yaml --link-files -r control.yaml

The comparison with the result obtained by refining the offsets initially
computed from the WCS information indicates that both methods lead to basically
the same result.

.. generada con --geometry 0,0,850,1200
.. convert combined_v3.png -trim combined_v3_trimmed.png
.. convert -delay 100 -loop 0 combined_v[13]_trimmed.png comparison_v3.gif


.. only:: html

   .. image:: comparison_v3.gif
      :width: 100%
      :alt: combined image, version 3 compared with version 1

.. only:: latex

   |pic_combined_v1_trimmed| |pic_combined_v3_trimmed|

   .. |pic_combined_v1_trimmed| image:: combined_v1_trimmed.png
      :width: 48%

   .. |pic_combined_v3_trimmed| image:: combined_v3_trimmed.png
      :width: 48%

.. generada con --geometry 0,0,850,1200 --bbox 1100,1600,800,1300
.. convert combined_v3_zoom.png -trim combined_v3_zoom_trimmed.png
.. convert -delay 100 -loop 0 combined_v[13]_zoom_trimmed.png comparison_v3_zoom.gif

.. only:: html
   .. image:: comparison_v3_zoom.gif
      :width: 100%
      :alt: combined image, version 3 compared with version 1

.. only:: latex

   |pic_combined_v1_zoom_trimmed| |pic_combined_v3_zoom_trimmed|

   .. |pic_combined_v1_zoom_trimmed| image:: combined_v1_zoom_trimmed.png
      :width: 48%

   .. |pic_combined_v3_zoom_trimmed| image:: combined_v3_zoom_trimmed.png
      :width: 48%

.. note::

   The conclusion of these comparisons is that the user can rely on the offsets
   computed from the WCS information in the image headers as a
   reasonable initial guess, but that these offsets need to be refined. Unless
   something is really wrong with that WCS information, the user probabily will
   not need to measure the offsets manually. Anyway, that option is always
   there just in case it is necessary.

.. _improving_skybackground_1:

Improving the sky background (problem #1)
-----------------------------------------

The first obvious way to improve the background computation is by masking the
objects present in the image. This masking process requires an initial object
detection, that must be carried out on the result of an initial combination.
For that reason, this masking requires to set ``iterations`` to a number larger
than zero. 

In addition, the user can indicate that the sky signal at each pixel must be
computed from the signal at the same pixel in a predefined number of images
(close in observing time).

The observation result file ``dithered_v4.yaml`` includes both options:

.. literalinclude:: dithered_v4.yaml
   :lines: 117-121
   :emphasize-lines: 2-3
   :linenos:
   :lineno-start: 117

Note that ``refine_offsets: True`` is also being used, but without setting
``offsets`` with an external ASCII file (i.e., the initial offsets will be
computed from the WCS information in the image headers).

Execute numina to start the reduction including object masking:

::

   (emir) $ numina run dithered_v4.yaml --link-files -r control.yaml

It is useful to subtract the new result from the one derived previously:

::

   (emir) $ numina-imath obsid_combined_v1_results/result_image.fits - \
      obsid_combined_v4_results/result_image.fits difference_v4.fits

.. generada con --geometry 0,0,850,1200
.. convert combined_v4.png -trim combined_v4_trimmed.png
.. convert difference_v4.png -trim difference_v4_trimmed.png
.. convert -delay 100 -loop 0 combined_v[14]_trimmed.png difference_v4_trimmed.png comparison_v4.gif

.. only:: html

   .. image:: comparison_v4.gif
      :width: 100%
      :alt: combined image, version 4 compared with version 1

.. only:: latex

   |pic_combined_v1_trimmed| |pic_combined_v4_trimmed|

   .. |pic_combined_v1_trimmed| image:: combined_v1_trimmed.png
      :width: 48%

   .. |pic_combined_v4_trimmed| image:: combined_v4_trimmed.png
      :width: 48%

   .. image:: difference_v4_trimmed.png
      :width: 48%
      :align: center

.. generada con --geometry 0,0,850,1200 --bbox 1100,1600,800,1300
.. convert combined_v4_zoom.png -trim combined_v4_zoom_trimmed.png
.. convert difference_v4_zoom.png -trim difference_v4_zoom_trimmed.png
.. convert -delay 100 -loop 0 combined_v[14]_zoom_trimmed.png difference_v4_zoom_trimmed.png comparison_v4_zoom.gif

.. only:: html

   .. image:: comparison_v4_zoom.gif
      :width: 100%
      :alt: combined image, version 4 compared with version 1

.. only:: latex

   |pic_combined_v1_zoom_trimmed| |pic_combined_v4_zoom_trimmed|

   .. |pic_combined_v1_zoom_trimmed| image:: combined_v1_zoom_trimmed.png
      :width: 48%

   .. |pic_combined_v4_zoom_trimmed| image:: combined_v4_zoom_trimmed.png
      :width: 48%

   .. image:: difference_v4_zoom_trimmed.png
      :width: 48%
      :align: center


.. _improving_skybackground_2:

Improving the sky background (problem #2)
-----------------------------------------

In all the previous examples, the combined images always exhibit variations in
the sky background that are clearly visible in the image borders. The reason
for that is that some individual exposures (in particular the first two
individual images), have a wrong image background. 

The problem is more severe in the regions where the number of images used for
the combination is lower:

.. convert combined_v4_npix.png -trim combined_v4_npix_trimmed.png
.. convert +smush 40 combined_v4_trimmed.png combined_v4_npix_trimmed.png data_npix_v4.png

.. image:: data_npix_v4.png
   :width: 100%
   :alt: combined data and number of pixels used in combination, version 4

.. only:: html

   .. image:: combined_v4_statistics.png
      :width: 100%
      :alt: combined data version 4, statistical analysis
   
.. only:: latex

   .. image:: combined_v4_statistics.png
      :width: 82%
      :alt: combined data version 4, statistical analysis
   
It is also possible to examine the sky-subtracted individual images (files
ending in ``_rfs_i?.fits`` within the ``work`` subdirectories):

.. cd obsid_combined_v4_work
.. numina-ximshow result_image_*_rfs_i1.fits --z1z2 [-200,300] --pdf skysub_v4.pdf --figuredict "{'figsize':(8, 10), 'dpi':100}"
.. convert -delay 50 -loop 0 skysub_v4.pdf skysub_v4.gif
.. convert skysub_v4.gif skysub_v4-%d.png

.. only:: html

   .. image:: skysub_v4.gif
      :width: 100%
      :alt: individual sky-subtracted images

.. only:: latex

   |pic_v4_0| |pic_v4_1| |pic_v4_2|
   |pic_v4_3| |pic_v4_4| |pic_v4_5|
   |pic_v4_6| |pic_v4_7| |pic_v4_8|
   |pic_v4_9| |pic_v4_10| |pic_v4_11|
   |pic_v4_12| |pic_v4_13|

   .. |pic_v4_0| image:: skysub_v4-0.png
      :width: 24%

   .. |pic_v4_1| image:: skysub_v4-1.png
      :width: 24%

   .. |pic_v4_2| image:: skysub_v4-2.png
      :width: 24%

   .. |pic_v4_3| image:: skysub_v4-3.png
      :width: 24%

   .. |pic_v4_4| image:: skysub_v4-4.png
      :width: 24%

   .. |pic_v4_5| image:: skysub_v4-5.png
      :width: 24%

   .. |pic_v4_6| image:: skysub_v4-6.png
      :width: 24%

   .. |pic_v4_7| image:: skysub_v4-7.png
      :width: 24%

   .. |pic_v4_8| image:: skysub_v4-8.png
      :width: 24%

   .. |pic_v4_9| image:: skysub_v4-9.png
      :width: 24%

   .. |pic_v4_10| image:: skysub_v4-10.png
      :width: 24%

   .. |pic_v4_11| image:: skysub_v4-11.png
      :width: 24%

   .. |pic_v4_12| image:: skysub_v4-12.png
      :width: 24%

   .. |pic_v4_13| image:: skysub_v4-13.png
      :width: 24%

One possibility is to remove the first two images from the list of images to be
reduced. This is undesirable because it obviously reduces the depth of the
combined image.

Another option is to apply an *ad hoc* correction, by fitting for example a
low-order 2D polynomial surface to the masked (i.e. removing objects)
sky-subtracted images. This option can be activated by using the requirement
``nside_adhoc_sky_correction``. We have incorporated that option in the
observation result file ``dithered_v5.yaml``, which also includes an iteration
to generate an object mask:

.. literalinclude:: dithered_v5.yaml
   :lines: 117-122
   :emphasize-lines: 5
   :linenos:
   :lineno-start: 117

Since the problem with the sky background in different for each quadrant of the
EMIR detector, the value of ``nside_adhoc_sky_correction`` indicates the
number subdivisions (in X and Y) in which each quadrant is subdivided. In this
case we are using a pattern of 10 x 10 regions in each quadrant. The median
value in each of these 100 subregions is computed (masking pixels affected by
objects) and a smooth spline surface is fitted to that collection of points.

::

   (emir) $ numina run dithered_v5.yaml --link-files -r control.yaml

.. generada con --geometry 0,0,850,1200
.. convert combined_v5.png -trim combined_v5_trimmed.png
.. convert -delay 50 -loop 0 combined_v[45]_trimmed.png comparison_v5.gif

.. only:: html

   .. image:: comparison_v5.gif
      :width: 100%
      :alt: combined image, version 5 compared with version 4

.. only:: latex

   |pic_combined_v4_trimmed| |pic_combined_v5_trimmed|

   .. |pic_combined_v4_trimmed| image:: combined_v4_trimmed.png
      :width: 48%

   .. |pic_combined_v5_trimmed| image:: combined_v5_trimmed.png
      :width: 48%

.. only:: html

   .. image:: combined_v5_statistics.png
      :width: 100%
      :alt: combined data version 5, statistical analysis

.. only:: latex

   .. image:: combined_v5_statistics.png
      :width: 82%
      :alt: combined data version 5, statistical analysis

In the last combination (v5) the sky background level is much flatter around
zero, except for those pixels in the combined image where only one single
exposure is available. By looking at the file ``result_i1_npix.fits`` it is
possible to check that those pixels are just at the borders of the combined
image.
