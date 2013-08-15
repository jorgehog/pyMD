from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

DCViz = "/home/jorgehog/code/DCViz/src"
numpy = "/usr/lib/python2.7/dist-packages/numpy/core/include"

setup(
    cmdclass = {'build_ext':build_ext},
    ext_modules = [Extension("pyMD_main", ["pyMD_main.pyx"],
                             include_dirs=[DCViz, numpy])]
)
