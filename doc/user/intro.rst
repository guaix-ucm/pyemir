
===================
PyEmir Installation
===================

This is PyEmir, the data reduction pipeline for EMIR. 

PyEmir is distributed under GNU GPL, either version 3 of the License, 
or (at your option) any later version. See the file COPYING for details.

PyEmir requires the following packages installed in order to
be able to be installed and work properly:

 - pyfits (http://www.stsci.edu/resources/software_hardware/pyfits)
 - setuptools (http://peak.telecommunity.com/DevCenter/setuptools)
 - numpy (http://numpy.scipy.org/)
 - scipy (http://www.scipy.org)

Additional packages are optionally required:
 - nose (http://somethingaboutorange.com/mrl/projects/nose) to run the tests
 - sphinx (http://sphinx.pocoo.org) to build the documentation

Webpage: https://guaix.fis.ucm.es/projects/emir

Maintainer: sergiopr@fis.ucm.es


Building 
--------

Run the following command::

   $ python setup.py build


For additional options::

   $ python setup.py build --help

Building the documentation
---------------------------
The Pyemir documentation is base on sphinx. With the package installed, the 
documentation can be built with a custom buld target::

  $ python setup.py build_sphinx


Installation
------------

Run the command::

   $ python setup.py install

For additional options::

   $ python setup.py install --help
