*************
Preliminaries
*************

.. warning::

   All the commands are assumed to be executed in a terminal running the **bash
   shell**.

   Don't forget to activate the same Python environment employed to install
   PyEmir.  In this document, the prompt ``(py36) $`` will indicate that this
   is the case.
   

Running PyEmir recipes from Numina
----------------------------------

The ``numina`` script is the interface with GTC pipelines. In order to execute
PyEmir recipes you should use execute something like:

::

   (py36) $ numina run <observation_result_file.yaml> -r <requirements_file.yaml>

where ``<observation_result_file.yaml>`` is an observation result in YAML
format, and ``<requirements_files.yaml>`` is a reqruirements file, also in YAML
format.

Note: YAML is a human-readable data serialization language (for details, see
`YAML Syntax desciption
<https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html>`_)

Use of interactive matplotlib plots
-----------------------------------

The interactive plots created by some Numina and PyEmir scripts have been
tested using the Qt5Agg backend of matplotlib. Nnote that this will require the
``pyqt`` package to be installed in your environment (see for example `the
conda instructions to install pyqt <https://anaconda.org/anaconda/pyqt>`_).

If you want to use the same backend, check that the following line appears in
the file ``.matplotlib/matplotlibrc`` (under your home directory):

::

   backend: Qt5Agg

If that file does not exist, generate it with the above line.

In most interactive matplotlib plots created by Numina and Pyemir you can press
``?`` over the plotting window to retrieve a quick help concerning the use of
some keystrokes to perform useful plot actions, like zooming, panning, setting
background and foreground levels, etc. Note that some of these actions are
already available in the navigation toolbar that appears at the top of the
plotting windows.


Installing ds9
--------------

Probably you already have ds9 installed in your system. If this is not the
case, you can use conda to do it!

::

   (py36) $ conda install ds9

Note that we have activated the ``py36`` environment prior to the installation
of the new package. That means that this particular ds9 installation will be
exclusively available from within that environment.


Installing dfits and fitsort
----------------------------

Two very useful utilities (used in this tutorial) are ``dfits`` and
``fitsort``. They will allow you to quickly examine the header content of a set
of FITS files. 

These two utilities belong to the ESO eclipse library. If you do
not have eclipse installed in your system, you can download the following
stand-alone files (somehow outdated, but they work perfectly fine for our
purposes and do not require anything else but a C compiler): 

- :download:`dfits.c <dfits.c>`
- :download:`fitsort.c <fitsort.c>`

These files can be directly compiled using any C compiler:

::

   $ cc -o dfits dfits.c
   $ cc -o fitsort fitsort.c

Note that it is highly advisable to place the two binary files in a directory
included in the path of your system.

