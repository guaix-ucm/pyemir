.. _mos_example:

***********
MOS example
***********

.. warning::

   All the commands are assumed to be executed in a terminal running the **bash
   shell**.

   Don't forget to activate the same Python environment employed to install
   PyEmir.  In this document, the prompt ``(emir) $`` will indicate that this
   is the case.

.. note::

   It is assumed that the reader has already followed the previous section of
   this tutorial :ref:`simple_example`. Some of the concepts already introduced
   there are not going to be repeated here with the same level of detail (or
   even mentioned at all!).

Let's consider the rectification and wavelength calibration of a MOS image with
slitlets configured in a non-longslit pattern. In this case, the scientific
observations are deep enough to allow the use of the airglow OH lines to
produce a good wavelength calibration. This means that no arc exposures are
needed!

Download the following file: `pyemir_mos_tutorial.tgz 
<http://nartex.fis.ucm.es/data/pyemir/pyemir_mos_tutorial_v1.tgz>`_.

Move to the directory where you have deployed the initial file tree structure
containing the basic PyEmir calibration files (see  :ref:`initial_file_tree`).

Decompress there the previously mentioned tgz file:

::

   (emir) $ tar zxvf pyemir_mos_tutorial_v1.tgz
   ...
   ...
   (emir) $ rm pyemir_mos_tutorial_v1.tgz

This action should have populated the file tree with 24 science exposures
(placed wihtin the ``data`` subdirectory) and some additional auxiliary files:

::

   (emir) $ tree
   ...

You can easily examine the header of the 12 science files using the utilities 
``dfits`` and ``fitsort`` (previously mentioned):

::

   (emir) $ dfits data/0002004* | fitsort object grism filter exptime date-obs
   FILE                                            	OBJECT       	GRISM   	FILTER  	EXPTIME   	DATE-OBS              	
   data/0002004226-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T00:40:52.36	
   data/0002004255-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T00:46:58.26	
   data/0002004284-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T00:53:02.05	
   data/0002004313-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T00:59:09.00	
   data/0002004342-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:05:30.72	
   data/0002004371-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:11:37.67	
   data/0002004400-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:17:41.46	
   data/0002004429-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:23:49.46	
   data/0002004458-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:30:12.23	
   data/0002004487-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:36:20.24	
   data/0002004516-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:42:24.03	
   data/0002004545-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:48:32.03	
   data/0002004574-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T01:54:54.80	
   data/0002004603-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:01:02.81	
   data/0002004632-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:07:06.59	
   data/0002004661-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:13:14.60	
   data/0002004690-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:19:37.37	
   data/0002004719-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:25:45.38	
   data/0002004748-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:31:49.16	
   data/0002004777-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:37:57.17	
   data/0002004806-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:44:19.93	
   data/0002004835-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:50:27.94	
   data/0002004864-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T02:56:31.73	
   data/0002004893-20190424-EMIR-STARE_SPECTRA.fits	goya_mask1a_J	J       	J       	359.986465	2019-04-25T03:02:40.79
 
It is also useful to display additional FITS keywords that store relevant
information concerning the observation strategy:

::

   (emir) $ dfits data/0002004*fits | fitsort  OBSBLOCK NIMGOBBL IMGOBBL NFRSEC FRSEC READMODE
   FILE                                            	OBSBLOCK	NIMGOBBL	IMGOBBL	NFRSEC	FRSEC	READMODE
   data/0002004226-20190424-EMIR-STARE_SPECTRA.fits	1       	4       	1      	28    	0    	RAMP    	
   data/0002004255-20190424-EMIR-STARE_SPECTRA.fits	1       	4       	2      	28    	0    	RAMP    	
   data/0002004284-20190424-EMIR-STARE_SPECTRA.fits	1       	4       	3      	28    	0    	RAMP    	
   data/0002004313-20190424-EMIR-STARE_SPECTRA.fits	1       	4       	4      	28    	0    	RAMP    	
   data/0002004342-20190424-EMIR-STARE_SPECTRA.fits	2       	4       	1      	28    	0    	RAMP    	
   data/0002004371-20190424-EMIR-STARE_SPECTRA.fits	2       	4       	2      	28    	0    	RAMP    	
   data/0002004400-20190424-EMIR-STARE_SPECTRA.fits	2       	4       	3      	28    	0    	RAMP    	
   data/0002004429-20190424-EMIR-STARE_SPECTRA.fits	2       	4       	4      	28    	0    	RAMP    	
   data/0002004458-20190424-EMIR-STARE_SPECTRA.fits	3       	4       	1      	28    	0    	RAMP    	
   data/0002004487-20190424-EMIR-STARE_SPECTRA.fits	3       	4       	2      	28    	0    	RAMP    	
   data/0002004516-20190424-EMIR-STARE_SPECTRA.fits	3       	4       	3      	28    	0    	RAMP    	
   data/0002004545-20190424-EMIR-STARE_SPECTRA.fits	3       	4       	4      	28    	0    	RAMP    	
   data/0002004574-20190424-EMIR-STARE_SPECTRA.fits	4       	4       	1      	28    	0    	RAMP    	
   data/0002004603-20190424-EMIR-STARE_SPECTRA.fits	4       	4       	2      	28    	0    	RAMP    	
   data/0002004632-20190424-EMIR-STARE_SPECTRA.fits	4       	4       	3      	28    	0    	RAMP    	
   data/0002004661-20190424-EMIR-STARE_SPECTRA.fits	4       	4       	4      	28    	0    	RAMP    	
   data/0002004690-20190424-EMIR-STARE_SPECTRA.fits	5       	4       	1      	28    	0    	RAMP    	
   data/0002004719-20190424-EMIR-STARE_SPECTRA.fits	5       	4       	2      	28    	0    	RAMP    	
   data/0002004748-20190424-EMIR-STARE_SPECTRA.fits	5       	4       	3      	28    	0    	RAMP    	
   data/0002004777-20190424-EMIR-STARE_SPECTRA.fits	5       	4       	4      	28    	0    	RAMP    	
   data/0002004806-20190424-EMIR-STARE_SPECTRA.fits	6       	4       	1      	28    	0    	RAMP    	
   data/0002004835-20190424-EMIR-STARE_SPECTRA.fits	6       	4       	2      	28    	0    	RAMP    	
   data/0002004864-20190424-EMIR-STARE_SPECTRA.fits	6       	4       	3      	28    	0    	RAMP    	
   data/0002004893-20190424-EMIR-STARE_SPECTRA.fits	6       	4       	4      	28    	0    	RAMP   
 
The last summary shows that:

- ``OBSBLOCK``: (counter of the block sequence) runs from 1 to 6

- ``NIMGOBBL=4`` (requested number of images per observing block): each block
  correspond to a typical ABBA dithering pattern along the slit.

- ``IMGOBBL`` (counter of the image sequence) runs from 1 to 4 within each
  block.

- ``NFRSEC=28`` (frames produced per detector exposure): note that
  ``READMODE=RAMP``, so 28 reads are performed in each exposure. Note that
  these frames are what GTC calls *raw* frames: each one of them corresponds to
  a single read in the ramp. These files are provided by GTC with the suffix
  ``*_raw.fits``. The images we are reducing here are the result of the fit of
  ramps to the collection of 28 reads (in each detector pixel). In this
  tutorial we are using the adjective *raw* to refer to the latter FITS images
  resulting from the ramp fitting process. 

- ``FRSEC=0`` (counter of the frame sequence): in this case, a value of ``0``
  indicates that the FITS image is the result of the ramp fitting process. In
  the individual FITS images ``*_raw.fits`` (not provided), corresponding to
  each detector read, this keyword runs from 1 to 28.

In summary, the 24 science images provided correspond to 6 ABBA blocks.
Have a look to any of these images.  For that purpose you can use ``ds9`` or
the visualization tool provided with numina:
   
::

   (emir) $ numina-ximshow data/0002004226-20190424-EMIR-STARE_SPECTRA.fits 

.. numina-ximshow data/0002004226-20190424-EMIR-STARE_SPECTRA.fits --geometry 0,0,650,850

.. image:: images/0002004226_raw.png
   :width: 800
   :alt: Image 0002004226 raw


The CSU configuration
=====================

It is clear from the previous figure that the EMIR slitlets were not configured
in `longslit` mode, but in `MOS` mode. In addition, it is important to
highlight that not all the slitlets were opened (i.e., the slits not assigned
to a particular scientific target were closed in order to avoid spurious
spectra in the image; note that even for the closed slitlets the corresponding
CSU bars are not completely closed to avoid collisions)

The slitlet configuration can be easily displayed with the help of the
auxiliary PyEmir script ``pyemir-display_slitlet_arrangement`` (please, note
the use in this case of the additional parameter ``--n_clusters 2``):

::

   (emir) $ pyemir-display_slitlet_arrangement \
     data/0002004226-20190424-EMIR-STARE_SPECTRA.fits \
     --n_clusters 2 --longslits
   ...
   ...
   ---> separator:   1.475

.. image:: images/0002004226_csu_configuration.png
   :width: 800
   :alt: Image 0002004226 csu configuration

The first figure displays the slitlet arrangement.

.. image:: images/0002004226_slitlet_widths.png
   :width: 800
   :alt: Image 0002004226 slitlet widths

Note that in this case we have employed the parameter ``--n_clusters 2``, that
forces the script to display a histogram of slitlet widths, and to compute two
clusters and a separating value (in this case 1.475 mm) which classifies the
slitlets in two groups: slitlets closed (widths < 1.475 mm) and opened (width >
1.475 mm). This number will be used later.


Preliminary rectification and wavelength calibration using OH lines
===================================================================

As explained in :ref:`simple_example`, the rectification and wavelength
calibration of any EMIR spectroscopic image can be obtained with two levels of
quality:

- **preliminary calibration**: without auxiliary calibration images, computed
  from the empirical calibration derived by the instrument team.

- **refined calibration**: that refines the empirical calibration by making use
  of either additional calibration images (i.e., arcs) or by using the airglow
  (OH) emission lines. The latter approach is going to be use in this example.
  
*The refinement process requires an initial estimation of the offsets in the
spatial (Y axis) and spectral (X axis) directions between the empirical
calibration and the actual data. These two offsets can be easily estimated
after computing the preliminary calibration.*

**At this point we are going to describe a new auxiliary script**
``pyemir-generate_yaml_for_abba`` **(not available in previous PyEmir versions)
which helps to create the observation result files required to execute the
reduction recipes. This is an excellent alternative to the manual edition and
writing of this kind of files.**

Before using the auxiliary script, one needs to generate a simple text file,
arbitrarily called ``list_abba.txt``,
with the list of scientific images (note that the list of files is generated
within the ``data`` subdirectory):

::

   (emir) $ cd data

   (emir) $ ls 0002004*fits > list_abba.txt

   (emir) $ cat list_abba.txt

   0002004226-20190424-EMIR-STARE_SPECTRA.fits
   0002004255-20190424-EMIR-STARE_SPECTRA.fits
   0002004284-20190424-EMIR-STARE_SPECTRA.fits
   0002004313-20190424-EMIR-STARE_SPECTRA.fits
   0002004342-20190424-EMIR-STARE_SPECTRA.fits
   0002004371-20190424-EMIR-STARE_SPECTRA.fits
   0002004400-20190424-EMIR-STARE_SPECTRA.fits
   0002004429-20190424-EMIR-STARE_SPECTRA.fits
   0002004458-20190424-EMIR-STARE_SPECTRA.fits
   0002004487-20190424-EMIR-STARE_SPECTRA.fits
   0002004516-20190424-EMIR-STARE_SPECTRA.fits
   0002004545-20190424-EMIR-STARE_SPECTRA.fits
   0002004574-20190424-EMIR-STARE_SPECTRA.fits
   0002004603-20190424-EMIR-STARE_SPECTRA.fits
   0002004632-20190424-EMIR-STARE_SPECTRA.fits
   0002004661-20190424-EMIR-STARE_SPECTRA.fits
   0002004690-20190424-EMIR-STARE_SPECTRA.fits
   0002004719-20190424-EMIR-STARE_SPECTRA.fits
   0002004748-20190424-EMIR-STARE_SPECTRA.fits
   0002004777-20190424-EMIR-STARE_SPECTRA.fits
   0002004806-20190424-EMIR-STARE_SPECTRA.fits
   0002004835-20190424-EMIR-STARE_SPECTRA.fits
   0002004864-20190424-EMIR-STARE_SPECTRA.fits
   0002004893-20190424-EMIR-STARE_SPECTRA.fits

   (emir) $ cd ..

We generate the observation result file to perform a preliminary rectification
and wavelength calibration:

::

   (emir) $ pyemir-generate_yaml_for_abba data/list_abba.txt \
     --step 0 --npreliminary 4 --outfile 0_rectwv_preliminary.yaml

   Expected sequence pattern: ABBA
   Number of sequences......: 6
   Full set of images.......: ABBAABBAABBAABBAABBAABBA
   --> File 0_rectwv_preliminary.yaml generated!

Note that the first argument after the name of the script in the text file
containing the list of scientific images. The rest of the arguments are:

- ``--step 0``: this number indicates the kind of observation result file that
  is going to be generated. The different options are
  ``0`` (initial rectification and wavelength calibration),
  ``1`` (refined rectificacion and wavelength calibration),
  ``2`` (fast reduction of ABBA observations), and
  ``3`` (careful reduction of ABBA observations).
  At this point we are interested in the initial calibration. The other options
  are going to be used later (see below).

- ``--npreliminary 4``: although we have a total of 24 scientific images 
  the initial rectification and wavelength calibration coefficients will be
  computed for the combined result of the first 4 images (i.e., the first ABBA
  sequence). We could use the 24 images in the combination but since it is
  possible that the images exhibit small offsets as a function of time, it is
  safer to use a smaller bunch of files. The OH lines in these images are 
  bright enough to provide a good wavelength calibration even in a single
  exposure.

- ``--outfile 0_rectwv_preliminary.yaml``: indicates the name chosen for
  the observation result file. This file (shown below) has a single block
  making use of the reduction recipe ``GENERATE_RECTWV_COEFF``:

  .. literalinclude:: 0_rectwv_preliminary.yaml
     :linenos:
     :lineno-start: 1

At this point, the preliminary rectification and wavelength calibration can be
carried out:

::

   (emir) numina run 0_rectwv_preliminary.yaml --link-files -r control.yaml
   ...
   ...


As expected, two new subdirectories have been created: 

- ``obsid_0002004226_rectwv_preliminary_results``

- ``obsid_0002004226_rectwv_preliminary_work``

::

   (emir) $ numina-ximshow \
     obsid_0002004226_rectwv_preliminary_results/reduced_mos.fits

.. numina-ximshow  obsid_0002004226_rectwv_preliminary_results/reduced_mos.fits --geometry 0,0,1060,850


.. image:: images/0002004226_preliminary.png
   :width: 800
   :alt: Image 0002004226 preliminary


Refined rectification and wavelength calibration
================================================

Alghough the rectified and wavelength calibration image that we have just
obtained appears to be fine, looking in detail it is possible to realize that
both the absolute and relative wavelength calibration between slitlets are
still not perfect, and that there exists a small offset between the expected
and the observed slitlet frontiers. Fortunately, both problems can be easily
solved.

Two approaches are possible:

1. **Manual determination of X and Y offsets:** We can proceed following
   the same process previously described in :ref:`simple_example`.

2. **Automatic determination of X and Y offsets:** The fact that the scientific
   images correspond to a real MOS observation, and not a longslit observation
   performed by alignining the slitlets, facilitates the use of 2D
   cross-correlation to estimate the offsets between the expected and the
   observed image distortions. Note that if the slitlets are aligned to
   simulate a longslit, the cross-correlation in the vertical direction is
   unable to determine hte vertical offset.

In this tutorial we are going to describe both approaches.

Manual determination of X and Y offsets
---------------------------------------

.. note::

   As described in :ref:`simple_example`, the task of visually finding the
   offsets can be performed with either the auxiliary PyEmir script
   ``pyemir-overplot_boundary_model``, or by using ``ds9`` with the auxiliary
   ds9-region files created during the preliminary rectification and wavelength
   calibration reduction. In the following two subsections we are using the
   latter option.

**Checking the spatial direction (Y axis)**

The offset in the spatial direction (Y axis) can be estimated by plotting the
expected slitlet frontiers (file ``ds9_frontiers_rawimage.reg`` placed within
the subdirectory ``obsid_0002004226_rectwv_preliminary_work``), derived in the
preliminary rectification and wavelength calibration, over the raw image:

::

   (emir) $ ds9 data/0002004226-20190424-EMIR-STARE_SPECTRA.fits &

- select ``scale --> zscale``

- select ``region --> load -->
  obsid_0002004226_rectwv_preliminary_work/ds9_frontiers_rawimage.reg``

.. image:: images/0002004226_overplot1.png
   :width: 800
   :alt: Image 0002004226 overplot 1

Zooming:

.. image:: images/0002004226_overplot2.png
   :width: 800
   :alt: Image 0002004226 overplot 2

From this visual examination one concludes that ``global_integer_offset_y_pix:
-3``. Note that the sign of the offset is chosen to place the actual data
within the predicted frontiers (displayed with dotted blue lines).

**Checking the wavelength direction (X axis)**

.. warning::

   The refinement process here described is based on the use of airglow (OH)
   emission lines in the science frames. This assumes that the target spectra
   are not dominant over the airglow emission. If this is not the case (for
   example when observing bright sources with short exposure times), the user
   should employ calibration arc images obtained before and/or after the
   science images. The refinement process should then be carried out as
   described in :ref:`simple_example`.

Continuing with the same ``ds9`` interface, overplot the expected location of
the airglow (OH) emission lines:

- select ``region --> load -->
  obsid_0002004226_rectwv_preliminary_work/ds9_oh_rawimage.reg``

.. image:: images/0002004226_overplot3.png
   :width: 800
   :alt: Image 0002004226 overplot 3

Note that the location of only the brightest OH lines are displayed. The visual
examination reveals that in this case ``global_integer_offset_x_pix: 1``. Note
that the sign of the offset is chosen to place the observed OH lines on the
predicted locations (displayed in cyan and magenta for the odd- and
even-numbered slitlets, respectively).

Improving the rectification and wavelength calibration
------------------------------------------------------

We are ready to obtain the refined rectification and wavelength calibration
coefficients using the integer offsets previously estimated. For that purpose
it is advisable to make use of the same auxiliary script that we used before:

::

   (emir) $ pyemir-generate_yaml_for_abba data/list_abba.txt \
     --step 1 \
     --refine_wavecalib_mode 12 \
     --minimum_slitlet_width_mm 1.475 \
     --maximum_slitlet_width_mm 3.0 \
     --global_integer_offset_x_pix 1 \
     --global_integer_offset_y_pix -3 \
     --rectwv_combined \
     --outfile 1_rectwv_combined.yaml

   Expected sequence pattern: ABBA
   Number of sequences......: 6
   Full set of images.......: ABBAABBAABBAABBAABBAABBA
   --> File 1_rectwv_combined.yaml generated!

Note that in this case the following arguments have been used:

- ``--step 1``: this indicates that the reduction recipe
  ``GENERATE_RECTWV_COEFF`` must carry out a refined estimation of the
  rectification and wavelength calibration coefficients.

- ``--refine_wavecalib_mode 12``: this indicates that the image corresponds to
  a science exposure, deep enough to detect OH sky lines, and that we are
  asking for a refinement of the wavelength calibration using that information.
  Note that if we were using an arc image, this parameter should be set to 2
  instead of 12 (as described in :ref:`simple_example`).

- ``--minimum_slitlet_width_mm 1.475`` and ``--maximum_slitlet_width_mm 3.0``:
  minimum and maximum slitlet widths (mm) for a slitlet to be considered as a
  scientific slitlet. Note that these numbers are compatible with the histogram
  of slitlet widths that we have obtained previously using
  ``pyemir-display_slitlet_arrangement`` with the parameter ``--n_clusters 2``.
  Only the slitlets which width was set within the specified range will be
  employed to derive a median sky spectrum (needed for the cross-correlation
  algorithm that is taking place during the refinement process). This allows to
  skip unused slitlets which could introduce noise in the cross-correlation
  method employed in the wavelength direction.

- ``--global_integer_offset_x_pix 1`` and ``--global_integer_offset_y_pix
  -3``: these are the offsets between the raw images and the expected
  empirical calibration, estimated as previously described.

- ``--rectwv_combined`` is just a flag (i.e. it is present or not, but has no
  associated value). **This flag is very important:** when present it indicates
  that all the images will be combined in a single image before computing the
  refined rectification and wavelength calibration coefficients.

- ``--outfile 1_rectwv_combined.yaml`` indicates the name chosen for the
  observation result file. This file (shown below) has a single block making
  use of the reduction recipe ``GENERATE_RECTWV_COEFF``:

.. literalinclude:: 1_rectwv_combined.yaml
     :linenos:
     :lineno-start: 1

This file contains a single block with the 12 images. The ``requirements``
section appears with the parameters employed when invoking the script 
``pyemir-generate_yaml_for_abba``.

Execute the reduction recipe:

::

   (emir) $ numina run 1_rectwv_combined.yaml --link-files -r control.yaml
   ...
   ...

As expected, two new subdirectories have been created:

- ``obsid_0002004226_rectwv_combined_work``

- ``obsid_0002004226_rectwv_combined_results``


We can have a look to the reduced image:

::

   (emir) $ numina-ximshow obsid_0002004226_rectwv_combined_results/reduced_mos.fits

.. numina-ximshow obsid_0002004226_rectwv_combined_results/reduced_mos.fits --geometry 0,0,1060,850

.. image:: images/0002004226_refined.png
   :width: 800
   :alt: Image 0002004226 refined

It is also possible to display the synthetic image, generated during the
execution of the reduction recipe, with the expected location of the airglow
(OH) lines (line intensities are normalized in the range from 0.0 to 1.0):

::

   (emir) $ numina-ximshow obsid_0002004226_rectwv_combined_work/expected_catalog_lines.fits \
     --z1z2 0,0.3

.. numina-ximshow obsid_0002004226_rectwv_combined_work/expected_catalog_lines.fits --geometry 0,0,1060,850 --z1z2 [-0.1,0.4]

.. image:: images/0002004226_expected_oh_lines.png
   :width: 800
   :alt: Image 0002004226 expected oh lines

In the ``obsid_0002004226_rectwv_combined_work`` subdirectory you can also find
a file named ``crosscorrelation.pdf`` which contains a graphical summary of the
cross-correlation process in the wavelength direction for each useful slitlet.
A few pages of that file are the following:

.. howto:
   convert -density 200 crosscorrelation.pdf[121] 0002004226_crosscorrelation0.png
   convert -density 200 crosscorrelation.pdf[122] 0002004226_crosscorrelation1.png
   convert -density 200 crosscorrelation.pdf[123] 0002004226_crosscorrelation2.png
   convert -density 200 crosscorrelation.pdf[151] 0002004226_crosscorrelation3.png


.. image:: images/0002004226_crosscorrelation0.png
   :width: 800
   :alt: cross-correlation example 0

.. image:: images/0002004226_crosscorrelation1.png
   :width: 800
   :alt: cross-correlation example 1

.. image:: images/0002004226_crosscorrelation2.png
   :width: 800
   :alt: cross-correlation example 2

.. image:: images/0002004226_crosscorrelation3.png
   :width: 800
   :alt: cross-correlation example 3


.. warning::

   The refined rectification and wavelength calibration has been saved in the
   file ``rectwv_coeff.json`` stored in the
   ``obsid_0002004226_rectwv_combined_results`` subdirectory. This file can be
   applied, as described below, to any raw image (with the same CSU
   configuration). Although the name of the file is the same in all the
   ``results`` subdirectories, there is no confusion because each one resides
   within a subdirectory with an unique name.

   The files ``rectwv_coeff.json`` are stored in JSON format, and contains all
   the relevant information necessary to carry out the rectification and
   wavelength calibration of any image.

   Note: JSON is an open-standard file format that uses human readable text to
   transmit data objects (see details in `JSON description
   <https://www.json.org/>`_).


Automatic determination of X and Y offsets
------------------------------------------

As previously mentioned, the integer offsets in X and Y that we have estimated
before with some effort, can be automatically determined by PyEmir when the
spectroscopic images correspond to MOS observations with the slitlets not
following a longslit pattern.

In order to use this option, copy the previous observation result file and
introduce the changes that are described next:

::

  (emir) $ cp 1_rectwv_combined.yaml 1_rectwv_combined_auto.yaml

  (emir) $ edit 1_rectwv_combined_auto.yaml
  ...

.. literalinclude:: 1_rectwv_combined_auto.yaml
     :linenos:
     :lineno-start: 1

Changes in ``1_rectwv_combined_auto.yaml`` in comparison with
``1_rectwv_combined.yaml``:

- Line 1: the suffix ``_auto`` has been added to the block ``id``

- Lines 33 and 34: both integer global offsets have been set to 0 (they are
  going to be estimated automatically!)

- Insert the new line 35, indicating ``global_integer_offsets_mode: auto``

Execute the reduction recipe:

::

   (emir) $ numina run 1_rectwv_combined_auto.yaml --link-files -r control.yaml
   ...
   ...

The reduction now includes the generation of a synthetic image with the
expected location of the OH lines. This image is stored in the
``obsid_0002004226_rectwv_combined_auto_work`` subdirectory with the name
``synthetic_raw_image.fits``.

.. howto:
   cd obsid_0002004226_rectwv_combined_auto_work
   numina-ximshow synthetic_raw_image.fits  --geometry 0,0,650,850 --z1z2 0,40
   cd ..

.. image:: images/0002004226_synthetic.png
   :width: 800
   :alt: synthetic images with OH lines

.. howto:
   cd obsid_0002004226_rectwv_combined_auto_work
   numina-ximshow synthetic_raw_image.fits  --geometry 0,0,650,850 --z1z2 0,40 --bbox 1100,1400,1100,1400 --> xxx1.png
   numina-ximshow reduced_image.fits  --geometry 0,0,650,850 --bbox 1100,1400,1100,1400 --z1z2 0,1200 --> xxx2.png
   convert -delay 100 -loop 0 xxx[12].png 0002004226_synthetic_comparison.gif
   cd ..
.. image:: images/0002004226_synthetic_comparison.gif
   :width: 800
   :alt: comparison between synthetic and original image

The automatic procedure has determined the following offset values:

- ``global_integer_offset_x_pix: 1``

- ``global_integer_offset_y_pix: -3``

These numbers are the same as the ones that we estimated previously. For that
reason, the resulting refined rectification and wavelength calibrations
coefficients are the same:

::

   (emir) $ diff \
     obsid_0002004226_rectwv_combined_auto_results/rectwv_coeff.json \
     obsid_0002004226_rectwv_combined_results/rectwv_coeff.json 

   8c8
   <   "uuid": "1139443e-99db-11e9-b393-3c15c2e3dc50",
   ---
   >   "uuid": "aca6d82e-99da-11e9-abf7-3c15c2e3dc50",
   18c18
   <     "creation_date": "2019-06-28T19:29:37.736426",
   ---
   >     "creation_date": "2019-06-28T19:26:49.004673",

Note that the contents of both files are identical except for the unique
``uuid`` (that identifies each particular calibration) and the
``creation_date``.


Fast ABBA combination
=====================

The reduction recipe `` ABBA_SPECTRA_FAST_RECTWV`` performs a fast ABBA
combination since it employs a single rectification and wavelength calibration
solution for all the individual images. For that reason, all the A images in
the ABBA sequences can be combined before rectification and wavelength
calibration. The same occurs with the B images. The resulting A minus B
combination is then rectified and wavelength calibrated.

**Note that this recipe assumes that the images in** ``data/list_abba.txt``
**follow an ABBA pattern**.

::

  (emir) $ pyemir-generate_yaml_for_abba data/list_abba.txt \
    --step 2 --rectwv_combined --outfile 2_abba_fast.yaml 

   Expected sequence pattern: ABBA
   Number of sequences......: 6
   Full set of images.......: ABBAABBAABBAABBAABBAABBA
   --> File 2_abba_fast.yaml generated!
   
The observation result file ``2_abba_fast.yaml`` has the following content:

.. literalinclude:: 2_abba_fast.yaml
     :linenos:
     :lineno-start: 1

The ``requirements`` section (lines 29-36) is now different:

- ``pattern``: indicates that the images in the ``frames`` sections correspond
  to ABBA sequences. That means that the total number of images must be an
  interger multiple of 4.

- ``rectwv_coeff``: the value that appears here is the expected file name for
  the rectification and wavelength calibration coefficients of the first image.
  Since the calibrations are assumed to be placed within the ``data``
  subdirectory, note that the path of this calibration starts with a `../` in
  order to indicate that the calibration is not in the default place but in the
  ``results`` subdirectory of that image.

- ``method``: combination method. The three available possibilities are
  ``sigmaclip``, ``mean`` and ``median``.

- ``method_kwargs``: additional parameters that control the behavior of the
  combination method. These options are only available for the ``sigmaclip``
  method, where the user can set the ``low`` and ``high`` rejection limits (in
  units of the standard deviation).

- ``voffset_pix``: indicates the offset (in pixels) between the A and B images.
  Leaving this number as zero means that no attempt to merge the A and B images
  is going to be performed.

::

  (emir) $ numina run 2_abba_fast.yaml --link-files -r control.yaml

Two results are produced by this recipe (see subdirectory
``obsid_abba_fast_results``):

- ``reduced_mos_abba.fits``: computed as the sum of all the A images minus the
  sum of all the B images. This image should exhibit A spectra with positive
  signal (white) and B spectra with negative signal (black).

- ``reduced_mos_abba_combined.fits``: in this image the coadded A and B spectra
  are combined by multiplying the coadded B spectra in
  ``reduced_mos_abba.fits`` by minus one (in order to obtain positiva signal),
  shifting the result along the spatial direction by ``voffset_pix`` pixels,
  and summing the coadded A spectra. When ``voffset_pix=0`` this image is
  empty.

.. numina-ximshow obsid_abba_fast_results/reduced_mos_abba.fits --geometry 0,0,650,850 --z1z2 [-5,5]

.. image:: images/0002004226_abba_fast.png
   :width: 800
   :alt: fast combination of ABBA images

Making a zoom in the resulting combined image, we estimate that the vertical
offsets between A and B spectra (white and black, respectively) is around 38
pixels. We can introduce this number in the requirement ``voffset_pix``:

::

   (emir) $ cp 2_abba_fast.yaml 2_abba_fast_bis.yaml

   (emir) $ edit 2_abba_fast_bis.yaml
   ...

We modify the label in line 1 (adding the suffix ``_bis``), and set the
vertical offset to 38 pixels in line 36:

.. literalinclude:: 2_abba_fast_bis.yaml
   :linenos:
   :emphasize-lines: 1, 36
   :lineno-start: 1

::

  (emir) $ numina run 2_abba_fast_bis.yaml --link-files -r control.yaml

.. numina-ximshow obsid_abba_fast_bis_results/reduced_mos_abba_combined.fits --geometry 0,0,650,850 --z1z2 [-5,5]

.. image:: images/0002004226_abba_fast_combined.png
   :width: 800
   :alt: fast combination of ABBA images (A on top of B)


Careful ABBA combination
========================

In the previous combination we have assumed that in all the A images the
targets were in exactly the same location in the slitlets. The same was assumed
for the B images.

However, it is possible to check whether this is actually the case by looking
in detail at some bright sources. For example, in these observations, 4 stars
were employed to align the CSU on sky. The brightest one was place in slitlet
number 2. Let's select all the A images and display a zoom around this slitlet
region:

::

   (emir) $ list_a=`dfits data/0002004*.fits | fitsort IMGOBBL | awk '$2 == 1 || $2 == 4 {print $1}'`

   (emir) $  numina-ximshow $list_a \
    --bbox 86,2047,9,70 \
    --z1z2 0,1100 \
    --pdffile brightspectrum_a.pdf \
    --keystitle ra,dec,rotang."RA: %s, DEC: %s, rotang: %5.2f"

The first command defines the list of images for which the keyword ``IMGOBBL``
is either 1 or 4 (A images). The second command generates a PDF file in which
a zoom around slitlet number 2 is displayed for all the selected images.

.. howto:
   list_a=`dfits data/0002004*.fits | fitsort IMGOBBL | awk '$2 == 1 || $2 == 4 {print $1}'`
   numina-ximshow $list_a \
    --bbox 86,2047,9,70 \
    --z1z2 0,1100 \
    --pdffile brightspectrum_a.pdf \
    --keystitle ra,dec,rotang."RA: %s, DEC: %s, rotang: %5.2f"
   convert -delay 100 -loop 0 brightspectrum_a.pdf brightspectrum_a.gif

.. image:: images/brightspectrum_a.gif
   :width: 800
   :alt: variation of brightspectrum of A images

It is clear that the actual location of the A spectra change along the
sequence. Unfortunately, there are no bright spectra in B images (the stars in
A images disappear when moving from A to B; this was done on purpose with the
idea o maximizing the number of slitlets devoted to scientific targets) to see
what is happening with the spectra in the B images (note that here we are
displaying the original FITS images, where the sky background is still
present). Anyway, it is reasonable to think that similar shifts should also
take place between different B images.

A careful combination method should take care of these offsets along the
spatial direction before combining the individual images.

PyEmir try to handle this problem with the reduction recipe
``ABBA_SPECTRA_RECTWV`` which in principle is able to:

- compute relative offsets between A images before their combination

- compute relative offsets between B images before their combination

- compute relative offsets between the combined A images and the combined B
  images before generating the final combined image

In addition, this recipe can use the same rectification and wavelength
calibration coefficients for all the individual images (as we have done
previously) or different calibrations for each one.

Including relative offsets between A and B images
-------------------------------------------------

Since, as we have just seen, it is possible that there are offsets between the
different A and B locations, a first thing one can fix is to align all the A
images before coadding them, and the same with the B images.

Without going into too many details, the recipe ``ABBA_SPECTRA_RECTWV``
computes an averaged spatial profile of a particular target (hopefully bright
enough) in the all the A images, and the same for an object in the B images.
Cross-correlating these spatial profiles it is possible to compute de relative
offsets of all the A images relative to the first A image. The same process is
carried out with the B images, using in this case the first B image as
reference. Note that bright OH sky lines are masked before computing the
spatial profiles. 

When the object selected for the A and B images is the same, the recipe is also
able to compute the relative offset between the combined A and the combined B
images.  However, sometimes the science targets are too faint to provide useful
spatial profiles for the cross-correlation method to work properly, and the
bright stars employed to align the slitlets on the sky are only present in the
A images (in order to maximize the number of slitlets allocated to science
targets). In the latter scenario, the relative offset between the combined A
images and the combined B images cannot be automatically computed.

A template of observation result file making using of this recipe can be
generated with the help of the auxiliary script
``pyemir-generate_yaml_for_abba`` with the argument ``--step 3``:

::

   (emir) $ pyemir-generate_yaml_for_abba data/list_abba.txt \
     --step 3 --rectwv_combined --outfile 3_abba_template.yaml 

   Expected sequence pattern: ABBA
   Number of sequences......: 6
   Full set of images.......: ABBAABBAABBAABBAABBAABBA
   --> File 3_abba_template.yaml generated!
   
The observation result generated, ``3_abba_template.yaml`` has the following content:

.. literalinclude:: 3_abba_template.yaml
   :linenos:
   :emphasize-lines: 1, 36-46
   :lineno-start: 1

The ``id`` label (line 1) is again different (to avoid overwriting the previous
reductions). The ``requirements`` section now contains a new subsection called
``refine_target_along_slitlet`` (lines 36-46), where specific parameters for
the computation of the relative offsets between images appear.  Some of them
are predefined, but others need to be inserted by the user (all the *To Be
Defined* ``TBD`` values):

- ``npix_removed_near_ohlines``: before computing the spatial profiles of the
  selected targets in the A and B images, the bright OH lines are masked. This
  parameter indicates the number of pixels (at each side of the expected
  location of the lines) that are masked.

- ``nwidth_medfilt``: after masking the bright OH lines, a median filter is
  applied in the spectral direction with a size given by this parameter.

- ``save_individual_images``: a value of 0 indicates that some intermediate
  images (e.g. the individual shifted images after computing the relative
  offsets with respect to the first A or B images) do not need to be saved in 
  the ``work`` subdirectory. This is useful to safe disk space. Setting this
  parameter to 1 allows those images to be saved.
  
- ``ab_different_target``: 

  - Setting this parameter to 0 means that the selected target to determine the
    offsets in the A and in the B images is the same. In that case the recipe
    will automatically compute the relative offsets between the A and B images
    using cross-correlation again.

  - Setting this parameter to 1 or -1 indicates that the target is different
    and the recipe will apply the offset computed from the WCS information in
    the image headers. The sign of this parameter indicates the relative
    position of A and B targets in the images: a value of 1 means that the B
    targets are on top of the A targets in the original images, whereas -1
    means that the B targets appears below the A targets.

- ``vpix_region_a_target``: pixel interval in the spatial direction enclosing
  the bright target in the A images. This should include not only the target
  itself but the sky region around it. **Note that these numbers must be
  measured in the rectified and wavelength calibrated images**.

- ``vpix_region_a_sky``: pixel interval in the spatial direction where the sky
  level can be estimated for the previous bright target. 

- ``vpix_region_b_target`` and ``vpix_region_b_sky`` play the same role s
  the previous two parameters, but for the bright B target.

- ``list_valid_wvregions_a``: list (python-style) of wavelength intervals to be
  employed to compute the average spatial profile of the bright A target. If
  this parameter does not appear in the observation result file, it means that
  the whole spectral range will be employed. This parameter is useful when the
  target is only bright in some wavelength ranges.

- ``list_valid_wvregions_b``: this parameter is analogous to the previous one,
  but for the bright B target.

In order to decide appropriate values for these parameters, it is useful to
examine the result of the previous reduction, overplotting the ds9 regions
corresponding to the boundaries of the different slitlets:

::

   (emir) $ ds9 obsid_abba_fast_bis_results/reduced_mos_abba.fits &

- select ``scale --> zscale``

- select ``region --> load --> obsid_abba_fast_bis_work/ds9_boundaries_rectified.reg``

::

   (emir) $ cp 3_abba_template.yaml 3_abba.yaml

   (emir) $ edit 3_abba.yaml
   ...


.. literalinclude:: 3_abba.yaml
   :linenos:
   :lines: 36-46
   :lineno-start: 36

In this example we do not have spectra of any of the bright 4 stars in the B
images. This means that we are going to choose the bright star in slitlet 2 for
the A images and the target in slitlet 24 for the B images. For that reason we
set ``ab_different_target: 1`` (note that we use ``1`` instead of ``-1``
because the B exposures show the spectra above the B exposures).

::

   (emir) $ numina run 3_abba.yaml --link-files -r control.yaml
   ...
   ...

The execution time of this recipe is larger than the time required by the
faster ``ABBA_SPECTRA_FAST_RECTWV`` recipe because now:

- every single image has to be rectified and wavelength calibrated
  independently, prior to the determination of the relative offsets between
  images

- after the offset determination (via cross-correlation) the images have to be
  shifted in the spatial direction (using fractions of pixels, which takes
  additional time)

.. numina-ximshow obsid_abba_results/reduced_mos_abba_combined.fits --geometry 0,0,650,850 --z1z2 [-5,5]

.. image:: images/0002004226_abba_combined.png
   :width: 800
   :alt: careful combination of ABBA images (A on top of B)



