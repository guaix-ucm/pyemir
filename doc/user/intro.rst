
*******************
PyEmir Installation
*******************

This is PyEmir, the data reduction pipeline for EMIR. 

PyEmir is distributed under GNU GPL, either version 3 of the License, 
or (at your option) any later version. See the file COPYING for details.

PyEmir requires the following packages installed in order to
be able to be installed and work properly:

 
 - `setuptools <http://peak.telecommunity.com/DevCenter/setuptools>`_
 - `numpy <http://numpy.scipy.org/>`_
 - `scipy <http://www.scipy.org>`_
 - `astropy <http://www.astropy.org>`_ >= 1.1
 - `matplotlib <http://matplotlib.org/>`_
 - `six <https://six.readthedocs.io//>`_
 - `numina <https://pypi.python.org/pypi/numina/>`_ >= 0.15
 - `photutils <https://photutils.readthedocs.io/en/stable/>`_
 - `sep <https://github.com/kbarbary/sep>`_
 - `scikit-image <http://scikit-image.org/>`_

Additional packages are optionally required:
 - `pytest <http://pytest.org>`_ to run the tests
 - `sphinx`_ to build the documentation

Webpage: https://guaix.fis.ucm.es/projects/emir

Maintainer: sergiopr@fis.ucm.es

Stable version
--------------

The latest stable version of PyEmir can be downloaded from  
https://pypi.python.org/pypi/pyemir

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

    $ git clone https://github.com/guaix-ucm/pyemir.git

And then installed following the standard procedure:::

    $ cd pyemir
    $ python setup.py install

Using conda
***********

Install and configure conda. Then install the dependencies (you can create an environment):::

    $ conda create --name emir python=3
    $ source activate emir
    $ (emir) conda install numpy scipy astropy matplotlib six scikit-image
    $ (emir) conda install -c astropy photutils
    $ (emir) conda install cython pyyaml

The latest development version of the emirdrp source code can be retrieved
using git. In addition, we will need the latest version of numina:::

    $ git clone https://github.com/guaix-ucm/numina.git
    $ git clone https://github.com/guaix-ucm/pyemir.git

Then, to build and install emirdrp:::

    $ (emir) cd numina
    $ (emir) python setup.py build
    $ (emir) python setup.py install
    $ (emir) cd ../emirdrp
    $ (emir) python setup.py build
    $ (emir) python setup.py install


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
