#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name = 'dreamBeam',
      version = '0.1',
      description = 'Measurement equation framework for interferometrY in Radio Astronomy.',
      author = 'Tobia D. Carozzi',
      author_email = 'tobia.carozzi@chalmers.se',
      packages = find_packages(), #['antpat', 'rime', 'telescopes'],
      license = 'BSD',
      classifiers = [
          'Development Status :: 1 - Planning',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: ISC License',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Visualization'
      ],
      install_requires=[
          'numpy>=1.10',
          'python-casacore',
          'matplotlib>=1.5',
          'antpat'
      ],
      scripts = ['scripts/pointing_jones.py']
     )
