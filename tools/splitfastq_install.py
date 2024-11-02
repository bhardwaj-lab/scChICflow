#!/usr/bin/env python3

from setuptools import setup
from Cython.Build import cythonize

setup(
    name='splitfastq',
    ext_modules=cythonize("splitfastq.pyx"),
)
