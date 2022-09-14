from setuptools import setup, Extension
from Cython.Build import cythonize


ext_modules = [
    Extension(
        "BinChillingTools",
        ['./src/BinChillingTools.pyx'],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
]

setup(
    name="BinChillingTools",
    ext_modules=cythonize(
        ext_modules,
        compiler_directives={'language_level' : "3", 'boundscheck': False})
)