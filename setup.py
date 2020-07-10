#!/usr/bin/env python

from setuptools import setup, find_packages


setup(name='pyemir',
      version='0.15.0',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='https://github.com/guaix-ucm/pyemir',
      license='GPLv3',
      description='EMIR Data Processing Pipeline',
      packages=find_packages(),
      package_data={
          'emirdrp.simulation': ['*.dat'],
          'emirdrp.instrument.configs': [
              'bars_nominal_positions_test.txt',
              'component-*.json',
              'instrument-*.json',
              'lines_argon_neon_xenon_empirical.dat',
              'lines_argon_neon_xenon_empirical_LR.dat',
              'Oliva_etal_2013.dat',
              'setup-*.json'

          ],
          'emirdrp': ['drp.yaml'],
      },
      test_suite="emirdrp.tests",
      install_requires=[
          'setuptools>=36.2.1', 'numpy', 'scipy',
          'numina>=0.22', 'astropy>=2',
          'enum34;python_version<"3.4"',
          'matplotlib', 'six', 'photutils>=0.2',
          'sep>0.5', 'scikit-image>=0.11', 'scikit-learn>=0.19', 'lmfit'
      ],
      extras_require={
          'tools':  ["PyQt5"],
          'docs': ["sphinx"],
          'test': ['pytest', 'pytest-remotedata']

      },
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
                  'emirdrp.tools.convert_refined_multislit_bound_param:main',
              'pyemir-display_slitlet_arrangement = ' +
                  'emirdrp.tools.display_slitlet_arrangement:main',
              'pyemir-fit_boundaries = ' +
                  'emirdrp.tools.fit_boundaries:main',
              'pyemir-generate_yaml_for_abba = ' +
                  'emirdrp.tools.generate_yaml_for_abba:main',
              'pyemir-generate_yaml_for_dithered_image = ' +
                  'emirdrp.tools.generate_yaml_for_dithered_image:main',
              'pyemir-median_slitlets_rectified = ' +
                  'emirdrp.processing.wavecal.median_slitlets_rectified:main',
              'pyemir-merge_bounddict_files = ' +
                  'emirdrp.tools.merge_bounddict_files:main',
              'pyemir-merge2images = ' +
                  'emirdrp.tools.merge2images:main',
              'pyemir-rectwv_coeff_add_longslit_model = ' +
                  'emirdrp.processing.wavecal.rectwv_coeff_add_longslit_model:main',
              'pyemir-rect_wpoly_for_mos = ' +
                  'emirdrp.tools.rect_wpoly_for_mos:main',
              'pyemir-rectwv_coeff_from_arc_image = ' +
                  'emirdrp.processing.wavecal.rectwv_coeff_from_arc_image:main',
              'pyemir-rectwv_coeff_from_mos_library = ' +
                  'emirdrp.processing.wavecal.rectwv_coeff_from_mos_library:main',
              'pyemir-rectwv_coeff_to_ds9 = ' +
                  'emirdrp.processing.wavecal.rectwv_coeff_to_ds9:main',
              'pyemir-select_unrectified_slitlets = ' +
                  'emirdrp.tools.select_unrectified_slitlets:main',
              'pyemir-slitlet_boundaries_from_continuum = ' +
                  'emirdrp.tools.slitlet_boundaries_from_continuum:main',
              'pyemir-overplot_boundary_model = ' +
                  'emirdrp.processing.wavecal.overplot_boundary_model:main',
              'pyemir-overplot_bounddict = ' +
                  'emirdrp.tools.overplot_bounddict:main',
              ],
      },
      classifiers=[
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3.5",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          'Development Status :: 3 - Alpha',
          "Environment :: Console",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License (GPL)",
          "Operating System :: OS Independent",
          "Topic :: Scientific/Engineering :: Astronomy",
      ],
      long_description=open('README.rst').read()
      )
