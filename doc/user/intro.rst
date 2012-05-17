
===================
PyEmir Installation
===================

This is PyEmir, the data reduction pipeline for EMIR. 

PyEmir is distributed under GNU GPL, either version 3 of the License, 
or (at your option) any later version. See the file COPYING for details.

PyEmir requires the following packages installed in order to
be able to be installed and work properly:

 
 - `setuptools <http://peak.telecommunity.com/DevCenter/setuptools>`_
 - `numpy <http://numpy.scipy.org/>`_
 - `scipy <http://www.scipy.org>`_
 - `pyfits <http://www.stsci.edu/resources/software_hardware/pyfits>`_
 - `pywcs <http://stsdas.stsci.edu/astrolib/pywcs/>`_
 - `matplotlib <http://matplotlib.sourceforge.net/>`_
 - `numdisplay <http://stsdas.stsci.edu/numdisplay/>`_

Additional packages are optionally required:
 - `nose <http://somethingaboutorange.com/mrl/projects/nose>`_ to run the tests
 - `sphinx`_ to build the documentation

Webpage: https://guaix.fis.ucm.es/projects/emir

Maintainer: sergiopr@fis.ucm.es

Stable version
--------------

The latest stable version of PyEmir can be downloaded from  
ftp://astrax.fis.ucm.es/pub/software/emir/

To install PyEmir, use the standard installation procedure:::

    $ tar zxvf pyemir-X.Y.Z.tar.gz
    $ cd pyemir-X.Y.Z
    $ python setup.py install
    
The `install` command provides options to change the target directory. By default
installation requires administrative privileges. The different installation options
can be checked with::: 

   $ python setup.py install --help
   
Development version
-------------------

The development version can be checked out with:::

    $ hg clone https://guaix.fis.ucm.es/hg/pyemir

And then installed following the standard procedure:::

    $ cd pyemir
    $ python setup.py install

Building the documentation
---------------------------
The PyEmir documentation is base on `sphinx`_. With the package installed, the 
html documentation can be built from the `doc` directory::

  $ cd doc
  $ make html
  
The documentation will be copied to a directory under `build/sphinx`.
  
The documentation can be built in different formats. The complete list will appear
if you type `make` 
  
.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _sphinx: http://sphinx.pocoo.org
