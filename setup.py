#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='razi',
      version='0.0.0',
      description='',
      #url='',
      packages = find_packages(exclude = [
            'docs',
            'docs.*',
            ]),
     )
