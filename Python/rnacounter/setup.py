

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

extensions = [
    Extension("rnacounter", ["rnacounter.pyx"],
              include_dirs=[np.get_include()],
             )
]

setup(
    name = "rnacounter",
    cmdclass = {'build_ext':build_ext},
    ext_modules = cythonize(extensions),
    #gdb_debug=True,
)

