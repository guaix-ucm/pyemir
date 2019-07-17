
.. _IMAGING_tutorial:

########################################################
Imaging mode tutorial: combination of dithered exposures
########################################################

.. warning::

   All the commands are assumed to be executed in a terminal running the **bash
   shell**.

   Don't forget to activate the same Python environment employed to install
   PyEmir.  In this document, the prompt ``(emir) $`` will indicate that this
   is the case.

This tutorial provides an easy introduction to the use of PyEmir (via Numina),
focusing on the combination of dithered exposures.

For detailed documentation concerning the installation of PyEmir, see the
:ref:`pyemir_installation` guide.

We strongly recommend to follow the different sections of this tutorial in the
provided order, starting with the simplest combination method, before
attempting to refine the combination procedure.

.. note::

   Two aspects that are still not covered by PyEmir in imaging mode are the 
   following:

   - Only integer offsets between images are considered: this is not especially
     important considering that the PSF is well oversampled. The benefit of
     using integer offsets is that the reduction speed is not increased.

   - Image distortions are not considered yet: this should be incorporated in a
     future release. Until them, the objects at the border of the combined
     images appear distorted (sorry!). The inclusion of this correction will
     imply the use of fractions of pixels for a proper image rectification
     prior to the combination, which will probably make the reduction more cpu
     demanding. Meanwhile, a temporary solution may be to combine the images
     using offsets between images computed using as references objects located
     in the image region where one requires the best image alignment.


.. only:: html

   **Tutorial index:**
   
.. toctree::
   :maxdepth: 2
   
   preliminary_combination
   refined_combination
   
