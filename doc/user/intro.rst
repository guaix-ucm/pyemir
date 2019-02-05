.. _pyemir_installation:

*******************
PyEmir Installation
*******************

PyEmir, the data reduction pipeline for EMIR, is a Python package for Python
2.7 and Python 3.5 or greater.

The easiest method of installing PyEmir is using prebuilt packages. You can
also build from the development version. 

Maintainers: sergiopr@fis.ucm.es, cardiel@ucm.es

.. warning::

   All the commands are assumed to be executed in a terminal running the **bash
   shell**


What method of installation should I use?
-----------------------------------------

- If you are familiar with one method, use it (conda or virtualenv), since both
  are fully supported.

- In macOS, there is a well-known compatibility problem between virtualenv and
  `matplotlib <https://matplotlib.org/faq/osx_framework.html>`_, so we recommend setting up conda.

- In Linux, virtualenv is easier to setup.

Install in virtualenv
---------------------

`Virtualenv <https:virtualenv.pypa.io/en/stable/installation/>`_ is a tool that
allows to create isolated Python environments.

Since Python version 3.3, there is also a module in the standard library called
``venv`` with roughly the same functionality.

The steps to install and run PyEmir within a virtual environment are:

1. **Create a virtual environment using either virtualenv or venv**

  In order to create a virtual environment called e.g. emir using ``venv``:

  ::
  
     bash$ python3 -m venv emir /path/to/emir

  With ``virtualenv``:
  
  ::

     bash$ virtualenv /path/to/emir

  The directory ``/path/to/emir`` represents the location of the environment.
  It can be any valid directory path.


2. **Activate the environment**

  After creating the environment, the directory ``/path/to/emir`` contains a
  Python tree. One of the directories is ``/path/to/emir/bin``, whih contains a
  script called ``activate``. To activate the environment we ``source`` (a bash
  shell command) this script file:

  ::
  
     bash$ source /path/to/emir/bin/activate

  which yields a different system prompt to the user:

  ::
  
     (emir) bash$

  Now, the name of the environment appears before the standard prompt. We can
  use the environment only on the consoles or terminals where we have
  previously activated it.

3. **Install PyEmir with pip**

  After the environment activation, we can install PyEmir with ``pip``. This is
  the standard Python tool for package management. It will download the package
  and its dependencies, unpack everything and compile when needed.

  ::
  
     (emir) bash$ pip install pyemir
     ...
     ...

4. **Test the installation**

  We can test the installation by running the ``numina`` command:

  ::

     (emir) bash$ numina
     DEBUG: Numina simple recipe runner version 0.19

5. **Update within the environment**

  In order to update PyEmir within a virtualenv installation the user should
  execute:
  
  ::
  
     (emir) bash$ pip install -U pyemir

6. **Deactivate the environment**
  
  To exit the environment is enough to exit the terminal or run the command
  ``deactivate``:

  ::
  
     (emir) bash$ deactivate
     bash$


Install in Conda
----------------

`Conda <https://conda.io/docs/>`_ was created with a target similar to
``virtualenv``, but extended its functionality to the management of packages in
different languages.

You can install `miniconda <https://conda.io/miniconda.html>`_ or `anaconda
<http://docs.anaconda.com/anaconda/install/>`_. The difference is that
miniconda provides a light-weight environment and anaconda comes with lots of
additional Python packages. By installing ``miniconda`` you reduce the amount
of preinstalled packages in your system (after installing ``miniconda`` it is
possible to install ``anaconda`` by executing ``conda install anaconda``).

If you have updated the ``$PATH`` system variable during the miniconda or conda
installation, you can call conda commands directly in the shell, like this:

::

   bash$ conda info

If not, you will need the add the path to the command, like:

::

  bash$ /path/to/conda/bin/conda info


In this guide we will write the commands without the full path, for simplicity.

Once conda is installed according to the corresponding miniconda or anaconda
instructions, the steps to execute and run PyEmir under conda are:

1. **Create a conda environment**

  With coda, environments are created in a centralised manner (under the
  subdirectory ``./envs`` in your conda tree), i.e., we do not pass the path to
  the environment:

  ::

     bash$ conda create --name emir

  The Pyhton interpreter used in this environment is the same version
  currently used by conda. You can select a different version with

  ::

     bash$ conda create --name emir python=3.6

2. **Activate the environment**

  Activate the environment:

  ::

     bash$ conda activate emir

  which yields a different system prompt to the user:

  ::

     (emir) bash$ 

3. **Install PyEmir with conda**

  After the environment activation, we can instal PyEmir using conda (we
  provide conda packages for PyEmir in the `conda-forge channel
  <https://conda-forge.org>`_):

  ::

     (emir) bash$ conda install -c conda-forge pyemir
     ...
     ...

4. **Test the installation**

  We can test the installation by running the ``numina`` command:

  ::

     (emir) bash$ numina
     DEBUG: Numina simple recipe runner version 0.19

5. **Update within the environment**

  In order to update PyEmir within the conda environment the user should
  execute:
  
  ::
  
     (emir) bash$ conda update pyemir

6. **Deactivate the environment**
  
  To exit the environment is enough to exit the terminal or run the following
  command:

  ::
  
     (emir) bash$ conda deactivate
     bash$


Installing the development version (using conda)
------------------------------------------------

The development version is the most updated working version of the code (use it
at your own risk!). For this version to work properly, some additinal python
packages must have been already installed in your system. 

In order to facilitate the installation of the additional packages, it is
useful to add the AstroConda channel:

::

   bash$ $ conda config --add channels http://ssb.stsci.edu/astroconda

It is easy to create a new environment and install the required
packages using (in this example python 3.6 is defined as the default python
interpreter):

::

   bash$ conda create --name emir python=3.7 \
   astropy \
   cython \
   ipython \
   jupyter \
   matplotlib \
   numpy \
   photutils \
   pytest \
   PyYaml \
   scikit-image \
   scipy \
   setuptools \
   six \
   sphinx

Activate the new environment:

::

   bash$ conda activate emir
   (emir) bash$

**Installing/updating numina**

Download and install the development version using git:

::

   (emir) bash$ git clone https://github.com/guaix-ucm/numina.git
   (emir) bash$ cd numina
   (emir) bash$ python setup.py build
   (emir) bash$ python setup.py install
   (emir) bash$ cd ..

If you have numina already installed in your system, but want to update the
code with the latest version, you need to move to the same directory where you
previously downloaded numina and reinstall it:

::

   (emir) bash$ cd numina
   (emir) bash$ git pull
   (emir) bash$ python setup.py build
   (emir) bash$ python setup.py install
   (emir) bash$ cd ..

Note: when updating numina, remember to update also pyemir (see next).

**Installing/updating pyemir**

After installing numina, you can install pyemir, following the same procedure
previously described for numina:

::
   
   (emir) bash$ git clone https://github.com/guaix-ucm/pyemir.git
   (emir) bash$ cd pyemir
   (emir) bash$ python setup.py build
   (emir) bash$ python setup.py install
   (emir) bash$ cd ..

If you have pyemir already installed in your system, but want to update the
code with the latest version, you need to move to the same directory where you
previously downloaded pyemir and reinstall it:

::

   (emir) bash$ cd pyemir
   (emir) bash$ git pull
   (emir) bash$ python setup.py build
   (emir) bash$ python setup.py install
   (emir) bash$ cd ..

Note: when updating pyemir, remember to update numina first (see above).

