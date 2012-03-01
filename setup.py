#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='pyemir',
      version='0.5.0',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/projects/emir',
      download_url='ftp://astrax.fis.ucm.es/pub/users/spr/emir/pyemir-0.5.0.tar.gz',
      license='GPLv3',
      description='EMIR Data Processing Pipeline',
      packages=find_packages('.'),
      package_data={'emir.simulation': ['*.dat'],
                     'emir.instrument': ['image_*.txt', 'spectrum_*.txt'],
                   },
      test_suite="nose.collector",
      install_requires=['setuptools', 'numpy', 'pyfits', 'scipy', 
		'sphinx', 'nose', 'pywcs',
		'matplotlib', 'numdisplay', 'numina>=0.5.0'],
      classifiers=[
                   "Programming Language :: Python",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
      long_description='''\
      This is PyEmir, the data reduction pipeline for EMIR.
      
      PyEmir is distributed under GNU GPL, either version 3 of the License, 
      or (at your option) any later version. See the file COPYING for details.
     
      PyEmir requires the following packages installed in order to
      be able to be installed and work properly:
      
      - setuptools (http://peak.telecommunity.com/DevCenter/setuptools)
      - numpy (http://numpy.scipy.org/)
      - scipy (http://www.scipy.org) 
      - pyfits (http://www.stsci.edu/resources/software_hardware/pyfits)
      - pyxdg (http://www.freedesktop.org/wiki/Software/pyxdg)
      - simplejson (http://undefined.org/python/#simplejson)
      
      Webpage: https://guaix.fis.ucm.es/projects/emir
      Maintainer: sergiopr@fis.ucm.es
            
      EMIR is a wide-field, near-infrared, multi-object spectrograph proposed 
      for the Nasmyth focus of GTC. It will allow observers to obtain from tens to 
      hundreds of intermediate resolution spectra simultaneously, in the 
      nIR bands Z, J, H and K. A multi-slit mask unit will be used for target acquisition. 
      EMIR is designed to address the science goals of the proposing team and 
      of the Spanish community at large. 
      ''',
      )
