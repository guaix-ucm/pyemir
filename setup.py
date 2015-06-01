#!/usr/bin/env python

from setuptools import setup, find_packages
import sys, os

# There is a problem installing/uninstalling with pip
# pip will uninstall pyemir AND numina 
# this is the bug https://github.com/pypa/pip/issues/355


sys.path.append(os.path.abspath('lib'))
import emir
version = emir.__version__
del sys.path[-1]

setup(name='pyemir',
      version=version,
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/projects/emir',
      license='GPLv3',
      description='EMIR Data Processing Pipeline',
      packages=find_packages('lib'),
      package_dir={'emir': 'lib/emir'},
      package_data={'emir.simulation': ['*.dat'],
                     'emir.instrument': ['image_*.txt', 'spectrum_*.txt'],
                     'emir': ['drp.yaml'],
                   },
      test_suite="emir.tests",
      install_requires=['setuptools', 'numpy', 'scipy',
                        'numina>=0.13.0', 'astropy>=0.4',
                        'matplotlib',
                        'scikit-image>=0.10'],
      use_2to3=True,
      zip_safe=False,
      entry_points = {
        'numina.pipeline.1': [
            'emir = emir.loader:load_drp',
            ],
        'numina.storage.1': [
            'emir_default = emir.loader:load_cli_storage',
            ]
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
