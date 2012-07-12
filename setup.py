#!/usr/bin/env python

from setuptools import setup, find_packages
import sys, os

BASE_PKGS=find_packages('src', exclude=['drp', 'drp.*'])
NAMESPACE_PKGS = ['numina.pipelines', 'numina.pipelines.emir']
ALL_PKGS = BASE_PKGS + NAMESPACE_PKGS

# There is a problem installing/uninstalling with pip
# pip will uninstall pyemir AND numina 
# this is the bug https://github.com/pypa/pip/issues/355


sys.path.append(os.path.abspath('src'))
import emir
version = emir.__version__

setup(name='pyemir',
      version=version,
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/projects/emir',
      download_url='http://astrax.fis.ucm.es/software/pyemir/pyemir-%s.tar.gz' % version,
      license='GPLv3',
      description='EMIR Data Processing Pipeline',
      packages=ALL_PKGS,
      package_dir={'emir': 'src/emir', 'numina.pipelines': 'src/drp'},
      package_data={'emir.simulation': ['*.dat'],
                     'emir.instrument': ['obsmodes.yaml', 
                     'image_*.txt', 'spectrum_*.txt'],
                   },
      test_suite="nose.collector",
      install_requires=['setuptools', 'numpy', 'pyfits', 
                        'scipy', 'sphinx', 'pywcs',
                        'matplotlib', 'numdisplay', 
                        'numina>=0.8.1'],
      test_requires=['nose',
                     ],
      # numdisplay lives here
      dependency_links = [
        'http://stsdas.stsci.edu/numdisplay'
        ],
      classifiers=[
                   "Programming Language :: Python :: 2.7",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
      long_description=open('README.txt').read()
      )
