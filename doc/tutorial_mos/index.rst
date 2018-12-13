
.. _MOS_tutorial:

###################
PyEmir MOS Tutorial
###################

This tutorial provides an easy introduction to the use of Numina and PyEmir,
focusing on the rectification and wavelength
calibration of EMIR spectroscopic images.

.. warning::

   This "MOS Tutorial" is still a work in progress; several aspects of PyEmir
   are not yet covered with sufficient detail. More instructions will be
   provided here in the future.

For detailed documentation concerning the installation of PyEmir, see the
:ref:`user`.

.. image:: pyemir_mos_scheme.png
   :width: 800
   :alt: PyEmir MOS Scheme

As shown in the previous diagram, PyEmir helps to generate a rectified and
wavelength calibrated 2D image. From this point, the astronomer can use her
favourite software tools to proceed with the spectra extraction and analysis.

.. image:: pyemir_mos_example.jpg
   :width: 800
   :alt: PyEmir MOS Scheme

The rectification and wavelength calibration of the original 2048x2048-size raw
images generates reduced 3400x2090-size images. The spatial and spectral
sampling of the raw images are preserved, as much as possible (to minimise
resampling problems) in the reduced images, and the spectral direction axis is
enlarged to acommodate the varying wavelength ranges covered by the slitlet
spectra depending on the location of the slitlets on the plane defined by the
CSU (Cold Slit Unit).

At present, PyEmir is able to work with raw images obtained with the following
spectroscopic configurations:

=====   ======
Grism   Filter
=====   ======
J       J
H       H
K       Ksp
LR      YJ
LR      HK
=====   ======

**Tutorial index:**

.. toctree::
   :maxdepth: 2

   preliminaries
   understanding
   simple_example
   mos_example
   ngc7798

