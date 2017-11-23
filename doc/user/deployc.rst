
============================
PyEmir Deployment with Conda
============================

Install conda
-------------

Install conda/anaconda from https://www.anaconda.com/download 
Both versions (2 and 3) are supported. 
If you have conda installed already, you don't need to do it again.

Create an environment
---------------------        
        
With the command::

    $ conda create --name emir python=3

Then, Activate the environment::

    $ source activate emir
    
Install dependencies
--------------------

Most of the dependencies can be grabbed from the conda repositories::

    $ conda install numpy scipy astropy matplotlib six scikit-image
    $ conda install cython pyyaml pytest
    $ conda install -c astropy photutils lmfit
    $ pip install sep
    
Download and install the source code
------------------------------------

The development version is hosted at Github.
Choose a top level directory for keeping pyemir source code, then::

    $ git clone https://github.com/guaix-ucm/numina.git
    $ git clone https://github.com/guaix-ucm/pyemir.git

Then, you can install the packges::

    $ cd numina
    $ python setup.py build && python setup.py install
    # lots of output
    $ cd ../pyemir
    $ python setup.py build && python setup.py install
    # lots of output
    $ cd ..

To check that the pipeline is installed, run::

    $ numina show-instruments

The expected output is::

    DEBUG: Numina simple recipe runner version 0.15.dev5
    Instrument: EMIR
     has configuration 'Default configuration' uuid=225fcaf2-7f6f-49cc-972a-70fd0aee8e96
    default is 'Default configuration'
     has datamodel 'emirdrp.datamodel.EmirDataModel'
     has pipeline 'default', version 1

