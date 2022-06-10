PyEmir Documentation
====================

Welcome. This is the Documentation for PyEmir (version |version|), 

EMIR_ is a **wide-field, near-infrared, multi-object spectrograph** (MOS)
installed at the Nasmyth focus of GTC_. Its MOS mode allows observers to obtain
**tens of intermediate-resolution spectra simultaneously** in the nIR bands *Y,
J, H* and *K*. EMIR is designed to address the science goals of the proposing
team and of the Spanish community at large.


Maintainers: Sergio Pascual (sergiopr@fis.ucm.es), and Nicolás Cardiel
(cardiel@ucm.es)


.. warning::

   If you are reducing EMIR data obtained with read-out mode RAMP, please,
   have a look to the instructions provided at `Checking the RAMPs
   <http://research.iac.es/proyecto/emir/pages/observing-with-emir/observing-utilities/checking-ramp-raw-data.php>`_.

   It is very easy to check the employed read-out mode using the auxiliary
   script ``fitsheader`` (provided in the astropy package). For example:

   ``$ fitsheader *.fits -k readmode -f``

.. only:: html

   Document index:

.. toctree::
   :maxdepth: 1
   
   installation/index
   Tutorials <https://guaix-ucm.github.io/pyemir-tutorials/index.html>
   user/index 
   reference/index
   glossary

.. _GTC: http://www.gtc.iac.es
.. _EMIR: http://www.gtc.iac.es/instruments/emir
.. _tut_prep: https://guaix-ucm.github.io/pyemir-tutorials/preliminaries/preliminaries.html
.. _tut_img: https://guaix-ucm.github.io/pyemir-tutorials/tutorial_imaging/index.html
.. _tut_spec: https://guaix-ucm.github.io/pyemir-tutorials/tutorial_mos/index.html
.. _tut_flat: https://guaix-ucm.github.io/pyemir-tutorials/tutorial_flat/index.html
