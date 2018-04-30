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
                        'numina>=0.16', 'astropy>=1.1',
                        'enum34;python_version<"3.4"',
                        'matplotlib', 'six', 'photutils>=0.2',
                        'sep>0.5', 'scikit-image>=0.11', 'lmfit'],
      zip_safe=False,
      entry_points={
          'numina.pipeline.1': [
              'EMIR = emirdrp.loader:load_drp',
          ],
          'console_scripts': [
              'pyemir-apply_rectification_only = ' +
                  'emirdrp.tools.apply_rectification_only:main',
              'pyemir-apply_rectwv_coeff = ' +
                  'emirdrp.processing.wavecal.apply_rectwv_coeff:main',
              'pyemir-continuum_flatfield = ' +
                  'emirdrp.tools.continuum_flatfield:main',
              'pyemir-convert_refined_multislit_param = ' +
                  'emirdrp.tools.convert_refined_multislit_param:main',
              'pyemir-display_slitlet_arrangement = ' +
                  'emirdrp.tools.display_slitlet_arrangement:main',
              'pyemir-fit_boundaries = ' +
                  'emirdrp.tools.fit_boundaries:main',
              'pyemir-merge_bounddict_files = ' +
                  'emirdrp.tools.merge_bounddict_files:main',
              'pyemir-rect_wpoly_coeff_add_longslit_model = ' +
                  'emirdrp.tools.rect_wpoly_coeff_add_longslit_model:main',
              'pyemir-rect_wpoly_for_mos = ' +
                  'emirdrp.tools.rect_wpoly_for_mos:main',
              'pyemir-rectwv_coeff_from_arc_image = ' +
                  'emirdrp.processing.wavecal.rectwv_coeff_from_arc_image:main',
              'pyemir-rectwv_coeff_from_mos_library = ' +
                  'emirdrp.processing.wavecal.rectwv_coeff_from_mos_library:main',
              'pyemir-select_unrectified_slitlets = ' +
                  'emirdrp.tools.select_unrectified_slitlets:main',
              'pyemir-slitlet_boundaries_from_continuum = ' +
                  'emirdrp.tools.slitlet_boundaries_from_continuum:main',
              'pyemir-overplot_boundary_model = ' +
                  'emirdrp.tools.overplot_boundary_model:main',
              'pyemir-overplot_bounddict = ' +
                  'emirdrp.tools.overplot_bounddict:main',
              'pyemir-verify_rect_wpoly = ' +
                  'emirdrp.tools.verify_rect_wpoly:main',
              ],
      },
      classifiers=[
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3.4",
                   "Programming Language :: Python :: 3.5",
                   "Programming Language :: Python :: 3.6",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
      long_description=open('README.rst').read()
      )
