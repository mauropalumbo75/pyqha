#!/usr/bin/env python

from distutils.core import setup, Extension

grunc = Extension('grunc', sources = ['grunc.c'])

setup(name = 'grunc',
      version = '1.0',
      description = 'Python Package with gruneisen C Extension',
      ext_modules = [grunc],
      author='Mauro Palumbo',
      author_email=''
     )
