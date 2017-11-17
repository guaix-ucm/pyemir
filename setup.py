#!/usr/bin/env python

from setuptools import setup, find_packages

# There is a problem installing/uninstalling with pip
# pip will uninstall pyemir AND numina 
# this is the bug https://github.com/pypa/pip/issues/355


import emirdrp
version = emirdrp.__version__

setup(name='pyemir',
      version=version,
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/projects/emir',
      license='GPLv3',
      description='EMIR Data Processing Pipeline',
      packages=find_packages(),
      package_data={
          'emirdrp.simulation': ['*.dat'],
          'emirdrp.instrument': ['image_*.txt', 'spectrum_*.txt'],
          'emirdrp.instrument.configs': [
              'instrument-225fcaf2-7f6f-49cc-972a-70fd0aee8e96.json',
              'component-1ac26dff-4c3f-43d8-a069-7470357120d6.json'
          ],
          'emirdrp': ['drp.yaml'],
      },
      test_suite="emirdrp.tests",
      install_requires=['setuptools', 'numpy', 'scipy',
                        'numina>=0.14', 'astropy>=1.1',
                        'matplotlib', 'six', 'photutils>=0.2',
                        'sep>0.5', 'scikit-image>=0.11', 'lmfit'],
      zip_safe=False,
      entry_points={
          'numina.pipeline.1': [
              'EMIR = emirdrp.loader:load_drp',
          ],
          'console_scripts': [
              'pyemir_apply_bpm = ' +
                  'tools.apply_bpm:main',
              'pyemir_fit_boundaries = ' +
                  'tools.fit_boundaries:main',
              'pyemir_display_slitlet_arrangement = ' +
                  'tools.display_slitlet_arrangement:main',
              'pyemir_merge_bounddict_files = ' +
                  'tools.merge_bounddict_files:main',
              'pyemir_slitlet_boundaries_from_continuum = ' +
                  'tools.slitlet_boundaries_from_continuum:main',
              'pyemir_overplot_boundary_model = ' +
                  'tools.overplot_boundary_model:main',
              'pyemir_overplot_bounddict = ' +
                  'tools.overplot_bounddict:main',
              ],
      },
      classifiers=[
                   "Programming Language :: Python :: 2.7",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
      long_description=open('README.rst').read()
      )
