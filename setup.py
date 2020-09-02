#!/usr/bin/env python

from setuptools import setup, find_packages
from dreambeam import __version__

setup(name='dreamBeam',
      version=__version__,
      description='Measurement equation framework for radio interferometry.',
      author='Tobia D. Carozzi',
      author_email='tobia.carozzi@chalmers.se',
      packages=find_packages(),
      package_data={'dreambeam.telescopes.LOFAR':
                    ['share/*.cc', 'share/simmos/*.cfg',
                     'share/alignment/*.txt'],
                    'dreambeam':
                    ['configs/*.txt']},
      include_package_data=True,
      license='ISC',
      classifiers=[
          'Development Status :: 1 - Planning',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: ISC License',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Visualization'
      ],
      install_requires=[
          'numpy>=1.10',
          'python-casacore',
          'matplotlib>2.0',
          'antpat>=0.4'
      ],
      entry_points={
        'console_scripts': [
            'pointing_jones = scripts.pointing_jones:cli_main',
            'FoV_jones = scripts.FoV_jones:main'
        ]

      }
      )
