#!/usr/bin/env python
from setuptools import setup, find_packages
from SiBreed.__init__ import __version__

setup(name='SiBreed',
      version=__version__,
      description='Breeding simulation program writen in Python.',
      author='Koki Chigira',
      author_email='chigirak@g.ecc.u-tokyo.ac.jp',
      url='https://github.com/KChigira/SiBreed/',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'matplotlib',
      ],
      entry_points={'console_scripts': [
            'parents = SiBreed.parents:main',
            'cross = SiBreed.cross:main',
            'selfing = SiBreed.selfing:main',
            'backcross = SiBreed.backcross:main',
            'randcross = SiBreed.randcross:main',
            'genostat = SiBreed.genostat:main',
            'genovisual = SiBreed.genovisual:main',
            'genomap = SiBreed.genomap:main',
            ]
      }
    )
