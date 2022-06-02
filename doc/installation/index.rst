.. _pyemir_installation:

*******************
PyEmir Installation
*******************

PyEmir, the data reduction pipeline for EMIR, is a Python package
(Python 3.6 or greater is required).

The easiest method of installing PyEmir is using prebuilt packages. You can
also build from the development version. 

Maintainers: Sergio Pascual (sergiopr@fis.ucm.es), and Nicolás Cardiel
(cardiel@ucm.es)

.. warning::

   Please, note that Windows is not a currently supported operative system.

.. warning::

   All the commands are assumed to be executed in a terminal running the **bash
   shell**.

Index:

- :ref:`pyemir_installation_what_method`

- :ref:`pyemir_installation_virtualenv`

- :ref:`pyemir_installation_conda`

- :ref:`pyemir_installation_development_version_venv`


.. _pyemir_installation_what_method:

What method of installation should I use?
-----------------------------------------

- If you are familiar with one method, use it (conda or virtualenv), since both
  are fully supported.

- In Linux, virtualenv is easier to setup (note: for still unknown reasons,
  some people using particular Ubuntu versions have faced some problems running
  --not installing-- the pipeline; we suggest those users to try the
  installation using conda).

- In macOS, there is a well-known compatibility problem between virtualenv and
  `matplotlib <https://matplotlib.org/faq/osx_framework.html>`_, so we
  recommend setting up conda.

**You do not need root privileges to install the software.
Everything will be installed under your home directory.**

We are explaining how to install the sofware using **environments**.
Environments allow the installation of software packages in isolated *work
places*, avoiding collisions between already installed software packages in
your computer. This approach has several advantages:

- You do not need to worry about using the correct python version because the
  environment will include the appropriate version.

- Additional python packages, required for the execution of the pipeline, will
  also be installed under your selected environment. If other versions of these
  packages are already installed in your computer, there will be no collisions
  between these different versions.

- Environments can be easily created and removed.

.. _pyemir_installation_virtualenv:

Install in virtualenv (or venv)
-------------------------------

`Virtualenv <https:virtualenv.pypa.io/en/stable/installation/>`_ is a tool that
allows to create isolated Python environments.

Since Python version 3.3, there is also a module in the standard library called
``venv`` with roughly the same functionality.

The steps to install and run PyEmir within a virtual environment are:

1. **Create a virtual environment using either virtualenv or venv**

  In order to create a virtual environment called e.g. emir using ``venv``:

  ::
  
     $ python3 -m venv emir /path/to/emir

  With ``virtualenv``:
  
  ::

     $ virtualenv /path/to/emir

  The directory ``/path/to/emir`` represents the location of the environment.
  It can be any valid directory path.


2. **Activate the environment**

  After creating the environment, the directory ``/path/to/emir`` contains a
  Python tree. One of the directories is ``/path/to/emir/bin``, whih contains a
  script called ``activate``. To activate the environment we ``source`` (a bash
  shell command) this script file:

  ::
  
     $ source /path/to/emir/bin/activate

  which yields a different system prompt to the user:

  ::
  
     (emir) $

  Now, the name of the environment appears before the standard prompt. We can
  use the environment only on the consoles or terminals where we have
  previously activated it.

3. **Install PyEmir with pip**

  After the environment activation, we can install PyEmir with ``pip``. This is
  the standard Python tool for package management. It will download the package
  and its dependencies, unpack everything and compile when needed.

  ::
  
     (emir) $ pip install pyemir
     ...
     ...

4. **Test the installation**

  We can test the installation by running the ``numina`` command:

  ::

     (emir) $ numina
     DEBUG: Numina simple recipe runner version 0.30

  The current PyEmir version can also be easily displayed using:

  ::

     (emir) $ numina show-instruments
     INFO: Numina simple recipe runner version 0.30.0
     Instrument: EMIR
      version is '0.16.0'
      has configuration 'Default configuration' uuid=225fcaf2-7f6f-49cc-972a-70fd0aee8e96
      has datamodel 'emirdrp.datamodel.EmirDataModel'
      has pipeline 'default', version 1
     

5. **Update within the environment**

  In order to update PyEmir within a virtualenv installation the user should
  execute:
  
  ::
  
     (emir) $ pip install -U pyemir

6. **Deactivate the environment**
  
  To exit the environment is enough to exit the terminal or run the command
  ``deactivate``:

  ::
  
     (emir) $ deactivate
     $

If at a given point you need to remove the environment, deactivate that
environment and delete the whole directory where the environment was created
(be careful with the use of this command; make sure you are deleting the
correct directory!):

::

   $ rm -fr /path/to/emir


.. _pyemir_installation_conda:

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

   $ conda info

If not, you will need the add the path to the command, like:

::

  $ /path/to/conda/bin/conda info


In this guide we will write the commands without the full path, for simplicity.

Once conda is installed according to the corresponding miniconda or anaconda
instructions, the steps to execute and run PyEmir under conda are:

1. **Create a conda environment**

  With coda, environments are created in a centralised manner (under the
  subdirectory ``./envs`` in your conda tree), i.e., we do not pass the path to
  the environment:

  ::

     $ conda create --name emir python=3

  Here we are asking that environment to be created including the last version
  of Python 3. If for any reason you need a particular Python version, you can
  specify it; for example, to force the use of Python 3.9:

  ::

     $ conda create --name emir python=3.9

2. **Activate the environment**

  Activate the environment:

  ::

     $ conda activate emir

  which yields a different system prompt to the user:

  ::

     (emir) $ 

3. **Install PyEmir with conda**

  After the environment activation, we can instal PyEmir using conda (we
  provide conda packages for PyEmir in the `conda-forge channel
  <https://conda-forge.org>`_):

  ::

     (emir) $ conda install -c conda-forge pyemir
     ...
     ...

4. **Test the installation**

  We can test the installation by running the ``numina`` command:

  ::

     (emir) $ numina
     DEBUG: Numina simple recipe runner version 0.30

5. **Update within the environment**

  In order to update PyEmir within the conda environment the user should
  execute:
  
  ::
  
     (emir) $ conda update pyemir

6. **Deactivate the environment**
  
  To exit the environment is enough to exit the terminal or run the following
  command:

  ::
  
     (emir) $ conda deactivate
     $

If at a given point you need to remove the environment, deactivate that
environment and remove it through conda:

::

   $ conda remove --name emir --all



.. _pyemir_installation_development_version_venv:

Installing the development version (using venv)
------------------------------------------------

The development version is the most updated working version of the code (use it
at your own risk!). 

::

   $ python3 -m venv venv_emir
   $ source venv_emir/bin/activate
   (venv_emir) $ git clone https://github.com/guaix-ucm/pyemir.git
   (venv_emir) $ cd pyemir
   (venv_emir) $ pip install -e .

