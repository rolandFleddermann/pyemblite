#!/usr/bin/env python

from setuptools import setup, find_packages

import numpy as np
from Cython.Build import cythonize

include_path = [np.get_include()]

ext_modules = cythonize(
        'pyemblite/*.pyx',
        include_path=include_path
    )
for ext in ext_modules:
    ext.include_dirs = include_path
    ext.libraries = ["embree3"]

setup(
    name="pyemblite",
    version='0.0.1',
    ext_modules=ext_modules,
    zip_safe=False,
    packages=find_packages(),
    package_data = {'pyemblite': ['*.pxd']}
)
