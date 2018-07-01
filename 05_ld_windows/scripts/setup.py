from distutils.core import setup
from Cython.Build import cythonize
import os

os.chdir("05_ld_windows/scripts/")
setup(ext_modules=cythonize("*.pyx"))
