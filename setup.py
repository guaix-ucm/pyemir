#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='pyemir',
      version='0.6.3',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/projects/emir',
      download_url='http://astrax.fis.ucm.es/software/pyemir/pyemir-0.6.3.tar.gz',
      license='GPLv3',
      description='EMIR Data Processing Pipeline',
      packages=find_packages('.'),
      package_data={'emir.simulation': ['*.dat'],
                     'emir.instrument': ['image_*.txt', 'spectrum_*.txt'],
                   },
      test_suite="nose.collector",
      install_requires=['setuptools', 'numpy', 'pyfits', 
                        'scipy', 'sphinx', 'pywcs',
                        'matplotlib', 'numdisplay', 
                        'numina>=0.6.0'],
      test_requires=['nose',
                     ],
      # numdisplay lives here
      dependency_links = [
        'http://stsdas.stsci.edu/numdisplay'
        ],
      data_files=[('share/numina/pipelines', ['emir-plugin.ini'])],
      classifiers=[
                   "Programming Language :: Python",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
      long_description=open('README.txt').read()
      )
