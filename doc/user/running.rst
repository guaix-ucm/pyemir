####################
Running the pipeline
####################

The EMIR DRP is run through a command line interface
provided by :program:`numina`.

The run mode of numina requires:
 
  * A observation result file in YAML_ format
  * A requirements file in YAML format 
  * The raw images obtained in the observing block
  * The calibrations required by the recipe
 
The observation result file and the requirements file are created by the user,
the format is described in the following sections.
 
********************************
Format of the observation result
********************************

The contents of the file is a serialized dictionary with the
following keys:

*id*: not required, string, defaults to 1
    Unique identifier of the observing block

*instrument*: required, string
    Name of the instrument, as it is returned by ``numina show-instruments``

*mode*: required, string
    Name of the observing mode, as returned by ``numina show-modes``

*frames*: required, list of strings
    List of images names

*children*: not required, list of integers, defaults to empty list
    Identifications of nested observing blocks

This is an example of the observation result file

.. code-block:: yaml

    id: dark-test-21
    instrument: EMIR
    mode: TEST6
    images:
      - r0121.fits
      - r0122.fits
      - r0123.fits
      - r0124.fits
      - r0125.fits
      - r0126.fits
      - r0127.fits
      - r0128.fits
      - r0129.fits
      - r0130.fits
      - r0131.fits
      - r0132.fits
   
*******************************
Format of the requirements file
*******************************

This file contains calibrations obtained by running recipes (called **products**)
and other parameters (numeric or otherwise) required by the recipes (named **requirements**). The file
is serialized using YAML_


Example requirements file:

.. code-block:: yaml

    version: 1                                             (1)
    products:                                              (2)
      EMIR:
       - {id: 1, content: 'file1.fits', type: 'MasterFlat', tags: {'filter': 'J'}, ob: 200}     (3)
       - {id: 4, content: 'file4.fits', type: 'MasterBias', tags: {'readmode': 'cds'}, ob: 400} (3)
      MEGARA:
       - {id: 1, content: 'file1.fits', type: 'MasterFiberFlat', tags: {'vph': 'LR-U'}, ob: 1200} (3)
       - {id: 2, content: 'file2.yml', type: 'TraceMap', tags: {'vph': 'LR2', 'readmode': 'fast'}, ob: 1203} (3)
    requirements: (4)
      MEGARA:
         default:
           MegaraArcImage:  (5)
              polynomial_degree: 5 (6)
              nlines: [5, 5]       (6)


1. Mandatory entry, ``version`` must be 1
2. Products of other recipes are list, by instrument
3. The products of the reduction recipes are listed. Each result must contain:
    * A ``type``, one of the types of the products of the DRP in string format
    * A ``tags`` field, used to select the correct calibration based on the keywords of
      the input.
    * A ``content`` field, a pointer to the serialized version of the calibration.
    * A ``id`` field, unique integer
    * A ``ob`` field, optional integer, used to store the observation id of the images that
      created the calibration.
4. Numerical parameters of the recipes are stored in ``requirements``, with different sections
   per instrument.
5. The name of the observing mode.
6. Different parameters for the recipe corresponding to the observing mode in (5)


********************
Running the pipeline 
********************

:program:`numina` copies the images (calibrations and raw data) from directory 
``datadir`` to directory ``workdir``, where the processing happens. 
The result is stored in directory ``resultsdir``. 
The default values are for each directory are ``data``, ``obsid<id_of_obs>_work`` and ``obsid<id_of_obs>_results``.
All these directories can be defined in the command line using flags::

  $ numina run --workdir /tmp/test1 --datadir /scrat/obs/run12222 obs.yaml -r requires.yaml

See :ref:`numina:cli` for a full description of the command line interface.

Following the example, we create a directory ``data`` in our current directory and copy
there the raw frames from ``r0121.fits`` to ``r0132.fits`` and the master bias ``master_bias-1.fits``.

The we run::

  $ numina run obsresult.yaml -r requirements.yaml
  INFO: Numina simple recipe runner version 0.15
  INFO: Loading observation result from 'obsrun.yaml'
  INFO: Identifier of the observation result: 1
  INFO: instrument name:
  ...
  numina.recipes.emir INFO stacking 4 images using median
  numina.recipes.emir INFO bias reduction ended
  INFO: result: BiasRecipeResult(qc=Product(type=QualityControlProduct(), dest='qc'), biasframe=Product(type=MasterBias(), dest='biasframe'))
  INFO: storing result

We get information of what's going on through logging messages. In the end, the result and log files are stored in ``obsid<id_of_obs>_results``.
The working directory ``obsid<id_of_obs>_work`` can be inspected too. Intermediate results will be saved here.


On the other hand, in the following we attach a short code to run pyemir
by using a Python script. This is useful to use the Python debugger.

.. code-block:: python

    from numina.user.cli import main

    def run_recipe():
        main(['run', 'obsresult.yaml', '-r', 'requirements.yaml'])

    if __name__ == "__main__":
        run_recipe()


.. _YAML: http://www.yaml.org
