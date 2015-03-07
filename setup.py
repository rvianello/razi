#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='razi',
      version='0.0.0',
      description='Using SQLAlchemy with chemical databases',
      classifiers=[
          "Development Status :: 1 - Planning",
          "Environment :: Plugins",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: MIT License",
          "Topic :: Scientific/Engineering :: Chemistry"
      ],
      url='http://razi.readthedocs.org/',
      keywords='chemistry cheminformatics sqlalchemy orm',
      packages=find_packages(exclude=[
          'docs',
          'docs.*',
      ]),
      install_requires=[
          'SQLAlchemy>=0.7.0',
      ],)
